## Native full rigid-frame evolution

`ForcedGrid3::try_new_with_rigid_frame_preserving_velocity` attaches the full
instantaneous frame force without an initial half kick. Serial/pooled
`advance_*_with_rigid_frame` and moving-wall variants evolve it in the original
collision/stream transaction. For cell-center displacement `r = x - origin`,
the non-Coriolis acceleration is
`b = g - a_origin - alpha cross r - Omega cross (Omega cross r)`.
The same admitted `RigidFrame3` position law and analytic skew inverse close
`F_next = rho_next*(b_next - 2*Omega_next cross u_next)` with
`rho_next*u_next = m_raw,next + F_next/2` after actual native transport and
wall reflection. Endpoint terms run in the existing tiles, without another
population buffer, source sweep or force iteration. The input acceleration
vector becomes accepted force storage only after full physical admission.

Caller-authored endpoint frame histories use lattice units and must supply
consistent angular velocity/acceleration and body-acceleration components in
frame axes. Grid and wall velocities are relative to those axes; no inertial
initial-condition remap is inferred. A periodic axis requires both Omega and
alpha parallel to it. Pure translation is compatible with all periodic axes;
zero frame motion retains the existing acceleration path exactly.

Original wall impulse, torque and relative wall work remain separate from
frame sources. Coriolis endpoint power vanishes; this does not assert a global
inertial energy balance or finite-inertia vessel dynamics. Refusal/cancellation
publishes no population, source history or step count. Focused regressions
exercise changing frames on D3Q19/D3Q27 collision families, actual transported
density/velocity, relative wall balances, free fall, zero-frame replay,
periodic-axis refusal, late Mach refusal and cancelled transaction replay.

## Transient reservoir work on the original source-caloric Stokes clock

The existing `march_natural_convection_with_enthalpy_and_viscous_heat`
entries, including implicit material endpoints and material solids, now admit
one opposed Pressure/Outlet pair alongside walls/symmetry. The physical-heating
opt-in selects the same normal-traction, zero-tangential-trace mouth operator
as steady pressure-work Stokes. Original spatial profiles, scalar programs,
independent plane histories and component-local pressure references remain
the mechanical pressure owners. Source thermodynamic pressure is unchanged.

The native backward-Euler momentum row remains intact. Each mirrored mouth
row represents a physical half cell: its kinetic storage and time loss use
one half of the mirrored `rho*area*width` inertia. Interior masses retain
their original arithmetic. Signed reservoir work is sampled independently
from the actual projected flux and the same centered boundary pressures used
by momentum. The transient sampler consumes already-centered solver pressure;
the steady public-field sampler subtracts its current physical datum exactly
once. Large changing common pressure cannot enter either work calculation.

`UnsteadyHeatingAudit::pressure_work_j` retains the signed contribution of
each original plane, and `supplied_pressure_work_j` retains the total.
Positive work supplies the fluid; negative work returns energy to a reservoir
during braking. Per-step and accumulated accounts close wall + body +
reservoir work + bounded continuity defect against physical dissipation,
actual native kinetic change and the separately reported backward-Euler loss.
Only original molecular/drag dissipation enters the existing caloric source.
Pressure supply, kinetic energy and numerical loss are never added as heat.

The same signed projected face field carries the original thermal inflow and
backflow enthalpy. Explicit source epochs, jointly converged implicit epochs,
nonlinear source h(T), material solids, contacts and thermostats retain their
existing owners. Failed admission, physical balance, source evaluation or
cancellation returns no accepted march or partial account. No new public
march wrapper, per-step pressure graph walk or pressure-field reconstruction
is introduced.

This is constant-density unsteady Stokes with the existing Reynolds ceiling
at most 0.1. Momentum convection and advected kinetic energy are omitted
together; the extension makes no Navier–Stokes open-energy claim. BDF2, SGS,
moving reservoir tangential traces, multiple opening axes, finite reservoir
storage and source EOS pressure evolution remain outside this heating owner.
The separate retained phase-momentum entry keeps its closed-domain admission.

Validation in this session: five tests compiled and passed using the exact
production kinetic-increment and extensive heating-account source, including
signed reservoir acceleration/braking and atomic balance/overflow refusal.
All eight changed Rust files parsed with the pinned formatter. Added native
regression definitions exercise graded-mouth inertia on all axes, changing
common-pressure invariance, and explicit/implicit source-caloric reversal with
physical heat and thermal backflow. These full native solver/thermal tests
were not compiled or executed because the fs-exec/asupersync dependency
closure is unavailable.

## Native implicit Coriolis evolution

`ForcedGrid3::try_new_with_coriolis_preserving_velocity` attaches the local
rotating-flow source without an initial half kick. Serial/pooled
`advance_*_with_coriolis` and moving-wall variants retain the accepted collision
force and resolve `F_next = rho_next*(a_next - 2*Omega cross u_next)` from the
actual native stream. The original analytic skew inverse solves physical
velocity and its half-force correction together. A constitutive readback check
refuses unresolved cancellation; native collision and final Mach limits remain.
The supplied acceleration vector becomes force storage in place, with no new
population buffer, force iteration or transport pass.

For fixed spin, isolated constant-density motion has Cayley angle
`-2*atan(|Omega|)` and zero Coriolis source work to roundoff; geostrophic balance
uses the same endpoint law. This is local rotating flow, without a frame-axis
remap or a transported fluid energy equation. Complete accelerating-frame
kinematics require the corresponding Euler/centrifugal/translational source.
Wall impulse, moment and work retain their native ownership. Source history,
populations and clock publish together after final admission/cancellation.
Zero spin retains existing acceleration arithmetic. Focused regressions cover
native D3Q19/D3Q27 families, inertial/geostrophic evolution, transported density,
moving-wall work, executor replay and refusal/retry.

## Independent pressure histories on native reservoir planes

`UnsteadyConfig::opening_pressure_programs` supplies an optional
`PressureOffsetProgram` on each original boundary plane in `Face3::ALL`
order. Each selected plane must be a Pressure/Outlet opening. Programs retain
one to sixteen finite signed-Pa knots inline, with strictly increasing
non-negative times, linear interpolation and endpoint holds. They add after
the existing global scalar pressure program or flow callback:
`p_face(t) = r0 + scale(t)*(p_authored - r0) + offset_side(t)`.
Thus different reservoirs may ramp and reverse independently while their
authored spatial profiles remain intact. Defaults preserve the existing path.

Each connected fluid component retains both its original datum `r0` and
the original opening plane that selected it. Equal original datum values do
not merge distinct anchor planes. At each endpoint the native coordinate
reference is `r0 + offset_anchor(t)`; the head used by momentum and reservoir
work is formed as
`scale*(p_authored-r0) + (offset_side-offset_anchor)`.
The offset difference is computed before adding the small authored head.
A large changing common pressure therefore cannot become acceleration or
erase a small head. Sealed components and face-fan ambient references keep
their original rules. Pressure readback and the physical-pressure continuity
work audit restore/subtract the same current component reference.

Both flow-only and shared coupled thermal/material loops install the complete
sampled load once before endpoint iteration. The same native solver, spatial
profiles, signed backflow fluxes, component topology and BE/BDF2 velocity
history persist through every trial. No per-step allocation or component
graph walk is introduced. Pressure is an instantaneous incompressibility
multiplier: changing its common datum shifts physical readback without
advancing momentum history or inserting a common pressure gradient.
Checkpointed momentum retains the old offset samples and rebases its pressure
guess into the receiving solver's current reference epoch, without restoring
a stale boundary load.

Admission checks all samples, doubled opening heads and physical references
before changing either the scalar factor or offsets. Invalid roles, overflow
or cancellation publish no partial load/march. Programs control at most six
whole boundary planes; separate openings on one plane share its program.
They prescribe incompressible mechanical pressure, not thermodynamic source
pressure, compressibility or reservoir storage. Selecting these controls alone
does not enable physical heating. The source-caloric viscous-heating opt-in
separately admits one opposed traction-mouth pair under its BE Stokes gates.

Focused regressions exercise production reference arithmetic, independent
all-axis flow reversal and thermal backflow under BE/BDF2, three-opening
equivalence with separately authored static profiles, disconnected components
with equal original datums but different anchor planes, checkpoint epochs,
native Stokes reservoir work/heat, large changing common pressure and atomic
refusal. The native integration definitions retain their original solver and
constitutive owners; no substitute flow or heat model is used.

Six bounded tests of the actual pressure-program and component-reference
arithmetic compiled and passed. The native solver, thermal/backflow,
checkpoint-owner and Stokes work/heat regression definitions were not compiled
or executed here because the complete dependency closure was unavailable.
Rust parsing and source review do not replace those integration gates.

## Acceleration endpoints during native sparse rigid-wall motion

`FreeSurface3::planar_rigid_wall_motion_with_body_acceleration_candidate` and
`advance_body_acceleration_with_planar_rigid_wall_motion` compose uniform/affine
body-acceleration endpoints with prescribed tangential translation and axial
spin of the existing completely flooded fixed plane. One native collision,
stream, gas reconstruction and conversion uses the accepted old force and
`F_next = rho_next * a_next(x)`. No zero-time reframe removes the changing load's
half impulse. Newly wet cells retain their own density and position. Actual
reflected populations still own wall impulse, angular impulse and mechanical
work; body acceleration is not added to those boundary receipts.

Every acceleration-owned native step now admits positive finite final wet
populations and actual force-corrected Mach after conversion, including cells
whose raw velocity alone would pass. All native/source/boundary/cancellation
refusals preserve the accepted state. Equal endpoints retain held-field bits.
Focused regressions cover heterogeneous transported density, endpoint impulse,
rigid-wall work, conversion, and refused high-Mach candidate replay. Motion is
prescribed over the interval; this fixed-plane isothermal path adds no normal
geometry motion, finite-inertia wall dynamics or fluid heat/energy equation.

## Retained pressure programs in coupled native marches

`UnsteadyConfig::pressure_program` now owns a validated inline
`PressureProgram` for native transient momentum **and** its shared coupled
thermal timeline. The default is exact factor-one identity. One to sixteen
finite `(time_s, factor)` knots are retained without heap allocation; time is
non-negative and strictly increasing. Factors may be positive, zero or
negative. Samples interpolate linearly between knots and hold the first/last
factor outside the tabulated interval. Each original step samples at its end
time; knots do not insert substeps. The total march end time must remain
finite.

The sampled factor acts on the original differential reservoir profile:
`p_face(t) = p_reference + factor(t) * (p_authored - p_reference)`.
Each connected fluid component keeps its original pressure reference and all
spatial face profiles. This drives the original momentum rows, pressure
projection and signed thermal backflow flux without rebuilding the solver,
walking the component graph or allocating during program evaluation. An
unchanged factor keeps the existing boundary-validation fast path. BE/BDF2
velocity history survives ramps and reversals; fan curves, body forces and
other loads retain their independent original inputs. A flow-only callback
from `simple_unsteady_with_pressure_schedule` or its SGS variant multiplies
the retained program and is still sampled exactly once per endpoint.

The shared application point precedes every thermal/material endpoint trial.
It covers ordinary and SGS conjugate marches, the wall-heat march, native
implicit Boussinesq feedback and the existing source-material explicit,
implicit, source-enthalpy and material-solid coupled entries that receive
`UnsteadyConfig`. Every nonlinear pass sees the same endpoint pressure load
and the same previous/older momentum levels. Original source identity,
thermodynamic pressure, composition, caloric law and heat storage ownership
are unchanged. Their existing validity and convergence refusals still apply.
Reference-density Boussinesq owners do not become compressible material flow.

The scalar program alone modulates one authored differential pressure profile;
independent per-plane histories additionally use `opening_pressure_programs`.
Neither evolves common-mode thermodynamic pressure. It does not add pressure control to a distinct liquid
owner that does not receive `UnsteadyConfig`, nor change physical-heating
admission: the current native transient viscous-heat entry still requires a
closed domain. Pressure forcing is not directly deposited as a heat source.

The native regression `tests/native_pressure_program.rs` checks actual
coupled SIMPLE/energy solves against independently assembled four-cell
upwind heat equations through pressure-driven reversal on all three axes,
including BE/BDF2 inertia, hot/cold reservoir exchange and closure. It also
checks large pressure-offset invariance, exact factor-one parity, composition
with the existing flow callback, retained source-density/viscosity feedback
with explicit and implicit material clocks, and refusal without partial
publication or input mutation. Module regressions cover knot admission,
endpoint holds and finite interpolation of extreme signed factors.

## Actual-density acceleration through native open reservoirs

`ForcedOpenBoundaryGrid3::try_new_with_body_acceleration_preserving_velocity`
attaches spatial lattice acceleration without changing physical velocity or
stress. Serial and pooled `advance_*_with_body_acceleration` accept a borrowed
next-epoch field. Collision retains `F_n`; each streamed donor and reconstructed
endpoint uses its own actual `F_(n+1)=rho_(n+1)*a_(n+1)`. At wall/open seams,
retained wall rows and the native half-source give a scalar affine density
closure, solved directly and checked against the final population mass.
Pressure/velocity profiles and D3Q19/D3Q27 source families retain their native meaning.

A reusable force vector keeps the authored acceleration immutable while other
cells still query their donors. The original two executor passes remain; each
stream tile resolves its endpoint force and physical Mach before publication.
The forced-open owner also exposes pooled held and evolving force-density steps.
Accepted populations, source history, step count and native wall receipt publish
together. Invalid fields, source closure, endpoint Mach or cancellation preserve
all accepted state. Body impulse is trapezoidal; reservoir mass/linear/angular
transport and prescribed-wall impulse/work remain distinct native contributions.

Focused regressions use heterogeneous density, all three reservoir axes,
profiled pressure/velocity, moving-wall seams/corners and native collision families;
independently converged prescribed-force steps check the direct closure.
Mass, linear/angular momentum, wall work, serial/pool identity, zero-source
bit preservation and refusal/retry are exercised. This is the existing fixed
geometry, low-Mach isothermal model; it adds no moving topology, implicit
velocity-dependent force law, impact acceleration or fluid energy equation.

## Scheduled differential pressure in native transient flow

`simple_unsteady_with_pressure_schedule` and its explicit SGS companion sample
one finite signed factor at each native time-step endpoint. Every authored
Pressure/Outlet value, including spatial profiles, becomes
`p_ref + factor * (p_authored - p_ref)`, using the existing fixed reference of
its own connected fluid component. Factor one preserves the original march;
zero removes reservoir head and negative values reverse it. This is scalar
modulation of the differential profile. Independent per-plane offsets use
`opening_pressure_programs`; neither control evolves a thermodynamic datum.

The original solver retains its BE/BDF2 velocity history, face inertia,
pressure projection, physical half-cell opening masses and backflow rules.
Internal fan curves, body loads and tangential reservoir velocities keep their
original owners. No per-step allocation or topology traversal is added:
changed samples only admit the original exposed pressure planes. Invalid or
unrepresentable samples and cancellation return no partial march and leave
caller-owned data unchanged. Discontinuous schedules require time steps aligned
with their authored transitions. This flow-only entry does not advance heat.

Four native regression definitions exercise all-axis pressure-driven inertia
and reversal with BE/BDF2, oscillatory temporal refinement and conservation,
spatial profiles with independent large component offsets and constant-one
parity, and rejected samples/cancellation. The changed Rust sources passed the
pinned compiler toolchain's syntax/format parser. Native execution remains
pending: the shared filesystem is full and no cached LBM test binary or complete
dependency outputs are available. No executed convergence claim is made here.

## Independent pressure datums in disconnected finite-volume flow

The native SIMPLE owner centers pressure independently in each face-connected
fluid component. Its first exposed pressure reservoir, in native face/plane
order, selects the numerical datum; a component reaching a face fan retains
that fan's ambient zero. Solids and unrelated components cannot choose its
datum. Pressure differences drive the original momentum and pressure-correction
equations. Reported cell pressure restores each component's authored reference,
while pressure-work heat uses the same local reference as momentum.

A large constant offset in a separate duct therefore cannot launch a spurious
startup flow, erase a small driving head, or overflow pressure admission across
regions that exchange no fluid. Connected/common-datum cases retain scalar
reference storage. Connectivity is resolved at initialization; no new
per-step graph traversal, pressure solve or matrix is introduced. Native
momentum checkpoints retain the reference map and rebase their pressure guess
when a changed phase mask joins or splits components.

Five focused regression definitions cover all-axis Poiseuille flow with
independent offsets, exact rest at unequal component pressures, pressure-work
and physical-heat invariance, checkpoint restoration/topology changes, and
fan-ambient ownership. Execution of these new regressions is pending because
the development executor became unavailable. These are incompressible
mechanical-pressure coordinates, not a gas EOS or a new pressure boundary law.
Result-changing conjugate semantics advance to version 10.

## Source-owned nonisothermal turbulent species mixing

`FvFlow::diffuse_turbulent_gas_species` connects the original momentum eddy
viscosity to finite `ClosedGasParcel` inventories with an explicit turbulent
Schmidt number. Canonical source species and typed mole-to-mass conversion
supply the fractions; the shared native graded-cell half-resistance supplies
`G`. All faces sample one original snapshot with
`dm_k = dt G (Y_left,k - Y_right,k)`. Exact absent species remain absent until
an actual incoming stream introduces their original source.

The existing material species-network owner transports each donor's partial
enthalpy, combines all streams, and reconstructs each affected fixed-volume
gas once. Its NASA caloric increments cancel arbitrary formation references
before thermal arithmetic. The returned source state supplies the updated
temperature, thermodynamic pressure, heat capacity, acoustic state, viscosity
and conductivity. No mixture-average enthalpy replacement, second EOS or
additional formation-energy heat is introduced.

Native cell IDs, dimensions, fluid mask, exact parcel volumes and accepted
carrier masses are admitted before stream construction. Common species on a
conducting face must retain compatible caloric and transport sources. The
original global withdrawal, local/global species and energy budgets remain
authoritative; source partitions, face zero-net-mass defects and aggregate
carrier-mass drift are also checked. Solid, exterior and zero-diffusivity
faces are impermeable. Refusal/cancellation publishes no partial endpoint.

This is a first-order explicit frozen-state component, retaining the original
256-parcel and 4096-stream limits. Excessive withdrawals refuse; incoming
material cannot replenish a donor within the same interval. Face conductances
and typed compositions are prepared once, and no per-face gas reconstruction
is performed. Time/mesh refinement remains necessary. Molecular diffusion,
turbulent heat flux, reactions and subsequent pressure/momentum evolution
retain their own operators. SIMPLEC gauge pressure is separate from the
returned thermodynamic pressure, and the carrier's old passive enthalpy is
not an alternate caloric state.

Focused regressions exercise independent unequal-temperature NASA integrals
on graded cells, formation-reference invariance, native Smagorinsky coupling,
Schmidt/heat-policy separation, introduction from exact species absence,
impermeable solids, exact null controls and immutable refusal.

## Spatial acceleration on shared D3Q19/D3Q27 transport

`ForcedGrid3::try_new_with_body_acceleration_preserving_velocity` initializes
`F = rho*a` in the native half-source frame without changing physical velocity
or stress. `advance_serial_with_body_acceleration` and its pooled/moving-wall
counterparts accept the next spatial acceleration in lattice units (`a_SI*dt²/dx`).
Collision uses the accepted force; each completed native stream tile resolves
`F_(n+1) = rho_(n+1)*a_(n+1)` using its actual transported, wall-corrected density.
Both the grid and native collision Mach limits apply to that physical endpoint.

The caller's vector becomes force storage in place. The original collision and
stream remain the only executor passes; workers perform endpoint admission and
retain canonical error ordering. Native wall receipts measure the same corrected
populations and retain their impulse/work meaning. Failure or cancellation
publishes neither populations, force history nor step count. Existing force-density
steps retain their meaning; a held acceleration requires explicit next endpoints.
The selected molecular/SGS collision field, fixed wall topology and low-Mach
isothermal scope remain native. These APIs do not add open or moving-geometry
boundary physics, implicit velocity-dependent loads, or a fluid energy equation.

## Plane-wall Mann turbulence at native D3Q19 inlets

With `turbulent-inflow`, `BoundaryMannInflow3::new_with_plane_wall` binds the
existing harmonic Mann slip-plane observation to the native D3Q19 velocity
profile. Point samples and exact face averages consume the complete corrected
field, including its evanescent wall term. Frozen snapshots, accepted OU
endpoints and exact accepted interval means retain the atmosphere owner's
modes, stochastic path and clock; applying a profile never advances the source.

The corrected field is divergence-free and has zero wall-normal velocity at
the authored plane. Mean flow and frozen translation must be tangent to that
plane. A mean-normal inlet-flux correction is admitted only when its direction
is also tangent. Every complete sampling footprint must remain above the
plane. Failed geometry, source matching, work/storage admission or cancellation
preserves the previously accepted profile.

Only one complete spatial filter is retained: point observations allocate no
gain array; finite faces store the original filter's required per-mode gains.
The same accepted profile feeds the original ordinary, pooled and SGS boundary
stepping. `spatial_filter()` now returns `Option<&MannBoxFilter>`; wall
observations instead expose `wall_spatial_filter()` for their complete
two-position covariance operator. A homogeneous sinc filter alone does not
represent wall-blocked covariance.

Focused regressions cover wall no-penetration and divergence through native
SGS stepping, independent spatial quadrature of finite faces, original OU
endpoint/interval clocks, and atomic refusal. This is outer-flow slip-plane
blocking, not a viscous no-slip layer or a wall-generated spectrum. Spatial
and temporal averages apply to velocity, not nonlinear momentum/energy fluxes.
The native low-Mach boundary and the authored turbulence model retain their
original resolution and applicability limits.

## Spatial body acceleration on native dense D3Q19

`BoundaryGrid3::rebind_body_acceleration_field` selects a prescribed local
acceleration in lattice units while retaining physical density, velocity and
force-corrected stress. `step_with_endpoint_body_acceleration_field` and its
pooled/selected-wall counterparts advance the original collision and stream
with `F_n`, then use each cell's transported density in `F_(n+1) = rho_(n+1)*a`.
The sampled local value follows that cell through donor observation, regularized
velocity/pressure reconstruction and the existing scalar wall-seam closure.
Body impulse is the trapezoidal sum of endpoint forces; actual wall/reservoir
transport remains separate.

Callbacks sample active cells once in x-fast order before mutation and skip
static solids. One acceleration vector becomes accepted force storage in place;
native population scratch is reused, with no extra collision or cloned grid.
Failed sampling, numerical admission or cancellation publishes neither source
nor populations. The prescribed field replaces both previous force components;
an ordinary step after rebinding holds force density and does not automatically
refresh acceleration. Velocity-dependent loads require their own implicit law.

G0/G1/G3/G4/G5 tests cover uniform-field parity, analytic transported-density
response, stress-preserving rebinding, local sampling, linear/angular wall
impulse, all six profiled reservoir orientations against independently
converged prescribed-force steps, serial/pooled identity, refusal and exact
retry. This extends the existing low-Mach lattice model; it does not establish
general compressible-flow validity or a whole-solver convergence certificate.

## Physical momentum at forced D3Q19 initialization

`BoundaryGrid3` constructors initialize unit-density fluid at physical rest in
the held Guo source frame, `f = f_eq(rho, 0) - S(0, F)/2`. Density perturbation
uses each cell's actual uniform-plus-spatial force. Thus initialization adds
neither a half-step body impulse nor the artificial `F F/(4 rho)` pressure
produced by shifting the equilibrium velocity. The native collision and stream
then supply their ordinary full-step momentum exchange. Boundary semantics
version 2 records this correction; zero-force population arithmetic and the
unforced boundary golden remain unchanged. `set_force_field` still explicitly
selects source coordinates; `rebind_force_field` preserves an evolved physical
state.

G0/G1/G3 regressions cover both force signs, all five D3Q19 collision families,
the initial isotropic pressure tensor, and integrated body impulse on periodic
grids with nonuniform density and spatial forcing. These establish the stated
initialization and balance laws, not a new whole-solver convergence claim.

The public tiled `Duct` constructor and its seeded density perturbation now
use the same physical-rest half-source convention. Both force signs and
nonuniform density therefore start with zero physical velocity; an interior
uniform-fluid probe gains one full force-density-over-density increment per
native step before wall signals arrive. Invalid perturbation amplitudes refuse
before mutation. Zero-force arithmetic, collision, streaming and wall rules
are unchanged. The historical core golden remains a kernel sentinel with its
original explicitly authored raw initial populations in a native unit test;
its hash is unchanged and does not certify the new initialization behavior.

## Positive Kelvin in physical conjugate heat transfer

The scalar and axis-aligned CHT entry points require positive finite Kelvin
for prescribed boundaries, ambient and backflow reservoirs, fixed cells, sinks
and initial storage states. Steady solves, transient steps and converged native
buoyancy results also admit their complete physical temperature fields and
compact junction temperatures before returning an accepted result. Cooling
that reaches zero or negative Kelvin refuses instead of publishing a clamped
or unphysical endpoint. CHT semantics version 9 records this admission.

Signed temperature increments and intermediate native k(T) linear proposals
remain numerical work. Their final physical endpoint must satisfy the Kelvin
domain; an intermediate proposal alone does not establish failure. Native
buoyancy pseudo-time and coupling candidates follow the same distinction.
The native lattice passive-scalar API retains its own scalar contract.
Existing normalized CHT benchmarks now declare an explicit positive Kelvin
datum and subtract it only from their analytic reference observations.

Focused regression definitions cover boundary and storage admission, steady
and transient overcooling, valid cooling with the original heat balance, and a
nonlinear conductive root whose first frozen-coefficient proposal is negative.

## Centered temperatures in conjugate conduction

The native scalar/axis-aligned finite-volume energy solve uses temperature
increments about a reference selected from held bath, prescribed-cell, sink or
old-time data. Boundary, sink and storage loads use temperature differences;
material conductivity, contacts, radiation and caloric queries keep physical
Kelvin. The original sparse operator and iteration budget remain in use.
Native constant-Cp advection retains its zero-K enthalpy convention through the
actual local volume-divergence term; source-owned affine enthalpy evaluates its
original tangent at the reference without replacing its caloric datum.

Physical face accounting centers the same face equation at its actual cell
temperature. Pure conduction therefore vanishes exactly at equal temperatures.
Nonlinear endpoint algebra checks subtract only the actual row-propagated
uncertainty from Kelvin rounding and coordinate subtraction, then apply the
requested increment residual tolerance. Independent physical heat checks remain
active; an absolute Kelvin-sized residual denominator cannot hide a small load.
These result-changing semantics are version 8 of `CHT_SEMANTICS_VERSION`.

Six focused regression definitions cover small heat loads across conductivity
contrasts and Kelvin datums, exact uniform equilibrium, finite-contact heat and
fixed-cell reactions, closed-body transient storage, and nonzero-divergence
open-flow accounting, including backflow-only reservoirs. Pinned Rust syntax
checks and five native tests of the original face kernels and reference selection
passed, including reciprocal affine-enthalpy transport. Full integration execution
remains pending after the dependency build exhausted shared resources; the native
kernel checks do not establish full solver, convergence or performance claims.

## Conducting walls in native open-flow heat transport

`ThermalOpenField3::set_wall_conduction(axis, links, controls)` attaches the existing implicit finite-wall graph to an original lower/upper external wall pair. Multiple axis pairs can coexist, including both transverse pairs of a rectangular duct. Link indices omit the selected normal coordinate and retain the original x-fast surface order. Every endpoint must face an active fluid cell; blocked wall-array slots cannot become hidden solid storage or consume sources.

Each graph prepares its reflections from actual post-collision D3Q7 populations before streaming. The original trial wall cells receive the accepted endpoints once; final joint publication retains the native flow/heat transaction. Rejected attempts alter only reusable graph scratch, which is rebuilt on retry. Internal conductive transfers enter the solid first-law account with opposite endpoint signs and do not become external heat. Prepared factors and workspace are reused; each pair has separate explicit graph resource bounds in addition to the carrier storage limit.

This enables forced throughflow with laterally conducting finite boundaries using one wall-temperature owner. Cross-axis graph connections, resolved solid depth, variable caloric properties, automatic buoyancy feedback and the full developing-channel correlation benchmark remain outside this implementation.

## Native open-flow temperature transport

`OpenBoundaryGrid3::prepare_thermal_field` binds a bounded D3Q7 scalar carrier to the original native D3Q19/D3Q27 geometry, solid mask, physical velocity and accepted time. `step_serial_with_thermal` stages the original momentum step and thermal transport with one final joint publication point. Sources use the existing equilibrium half-kicks and retained advective-flux history; optional per-cell thermal relaxation is admitted for the step. Persistent candidate arrays and histories are reused.

Every velocity/pressure face requires explicit reservoir temperature, including nominal outlets. Accepted local force-corrected normal velocity selects Danckwerts total-flux inflow or zero-conductive-flux outflow independently at each boundary link. Receipts measure actual incoming-minus-outgoing scalar populations. Stationary external walls retain independent finite cells per active fluid link; unassigned internal obstacle laws remain stateless and do not create hidden solid storage. Temperature, velocity, source, resource, cancellation and independently measured fluid/solid balance gates precede publication of populations, wall state, histories, receipts and counters.

The conserved scalar uses fixed reference volumetric heat capacity. This is passive low-Mach temperature transport, not compressible mass-specific enthalpy or automatic buoyancy feedback. Its heat receipt must not be replaced by the momentum reservoir mass receipt times temperature. The scalar component velocity envelope is additional to native Mach admission. No nonreflecting-outlet or general coupled time-order claim is made.

### Prescribed spatial and evolving force in native convection

`ForcedOpenBoundaryGrid3::prepare_thermal_field`, `step_serial_with_thermal`
and `advance_serial_with_thermal` carry the same original D3Q7 thermal owner
alongside native spatial force epochs. Initial thermal equilibrium and heat
advection use each cell's accepted physical velocity, including its local
half-force correction. Collision retains the opening force; reservoir
reconstruction and the final source half-kick use the actual endpoint force.
The tensor-diffusion variants retain the original full SPD diffusion kernel.

The final publication advances momentum, force, temperature, finite walls,
flux history, exterior-wall momentum/work receipts and counters together.
An accepted changing-force step exchanges the caller's endpoint buffer with
the old force buffer for reuse; refusal or cancellation preserves accepted
states, the last wall receipt and the caller's buffer. There is no per-step
adapter allocation, extra force reframe, inferred body-work heat or duplicate
thermal source. Reference heat capacity and the original opening-velocity
thermal split remain unchanged; no new coupled time-order claim is made.
Focused G1/G3 regressions cover analytic affine acceleration and source moments
on both native lattices, spatial heat-packet displacement, independent combined
fluid/wall heat closure, late refusal with exact retry, equivalent tensor and
held-force representations, and moving-wall receipt publication/cancellation.

## Finite immersed bodies in native convection

`ThermalOpenField3::set_immersed_thermal_bodies` groups explicitly selected
original solid voxels into finite lumped thermal bodies. Each body has ONE
temperature, heat capacity and held external source across every exposed
D3Q7 link. The original native flow mask and no-slip obstacle stay unchanged.
The same capability is available through `BoussinesqOpenGrid3`, retaining
native buoyancy/heat publication and the separate fluid density envelope.

Each body's backward-Euler balance sums the original reflected-pair resistance
law over actual post-collision outgoing populations. All faces use the common
endpoint; neither capacity nor source is multiplied by voxel/face count.
Per-body receipts measure actual reflected heat and independently computed
storage. Internal-wall and combined fluid/solid receipts include each body
once. Refusal or cancellation preserves body temperatures, accepted receipts,
thermal populations/history and native flow. Interface lists are prepared in
streaming order, with explicit additional storage bounds and no step allocation
or geometry lookup.

The independent-body installer uses an isothermal-interior approximation with fixed heat capacity and a
common authored bulk-to-surface resistance per exposed link. Its applicability
requires negligible internal temperature gradients; it does not resolve solid
diffusion, body contact, moving interfaces, variable caloric properties or
surface-area error of stair-stepped native voxels. Focused regressions check
the analytic ten-face/two-voxel step, capacity-weighted closed-box equilibrium,
separate hot/cold body inventories, periodic seams, exact refused-step retry,
and ownership, storage-budget and temperature-domain admission.

## Spatial immersed-solid conduction in native convection

`set_immersed_thermal_network` installs original immersed thermal nodes and
authored reciprocal conductances together. Different nodes retain distinct
temperatures, capacities and sources; a fully buried node may have zero fluid
interfaces. Each positive-conductance component must reach an exposed node.
The same installation is available through `BoussinesqOpenGrid3` without
detaching its force history or thermal publication.

One simultaneous backward-Euler solve couples all node storage, external
sources, solid conduction and actual post-collision fluid traces. The existing
`fs-sparse` LDLT implementation factors the diagonally scaled operator once;
later steps reuse the factor and retained linear workspace. No temperature
copy becomes a separate physical owner. Each graph edge supplies the same
endpoint heat with opposite signs to its two original nodes. Actual reflected
fluid heat, original capacity times representable temperature change, source
and conductive heat must close per node and over the complete graph. The
existing final thermal/flow publication is the only acceptance point.

Node/interface and graph/factor storage have separate explicit budgets.
Repeated undirected paths must be combined before installation; invalid
graphs, exhausted factor budgets, unsupported temperatures, numerical balance
failure and cancellation preserve the previous accepted state and receipts.
Geometry and fixed material/contact conductances remain caller-owned. This is
spatial lumped-cell conduction with first-order solid time accuracy; it does
not introduce nonlinear caloric laws, moving geometry, an automatic accuracy
certificate or a new stability guarantee for the fluid scheme. Focused tests
cover independent two-node balances including stiff coupling, heat release
from a fully buried core through an exposed shell, whole-capacity equilibrium,
late cancellation/exact retry and graph admission.

### Cartesian material voxels on the immersed network

`set_immersed_thermal_voxels` prepares the same native thermal network directly
from selected original solid cells. Each voxel retains its own fixed capacity,
temperature, isotropic conductivity and external source. Neighboring selected
voxels share the harmonic conductance of their two half cells; the fluid-facing
resistance adds the solid half-cell path to an explicitly authored contact
resistance. Interior solid contacts are ideal. No temperature or heat inventory
is maintained by the returned cell-to-node mapping.

The installer uses canonical native cell order, the original periodic axes and
actual Cartesian faces. Distinct periodic faces between the same two nodes add
in parallel; self-neighbors carry no heat. Unselected solid neighbors and
nonperiodic exterior solid faces insulate. The geometry preparation has its own
payload bound in addition to retained node/interface and sparse-factor bounds.
All installation checks precede replacement of the accepted network. The same
entry point is exposed through `BoussinesqOpenGrid3`.

Four executed native tests cover the independent heterogeneous two-cell balance,
half-solid surface resistance, periodic parallel faces and input-order invariance,
a sourced 3x3x3 solid with a buried center, and transactional topology/budget
refusal. Together with the graph, body and original owner tests, 38 native
tests passed against the actual sparse solver. Momentum-wrapper and complete
Cargo integration remain separate verification obligations. Spatial/time
refinement, nonlinear caloric laws, rotated anisotropy and subcell surface
accuracy are not established by these local balances.

The `immersed_body_cooling` example now accepts `resolved` (default) or `lumped`
solid treatment on the same native flow geometry. The resolved case applies
heat to the interior voxel and reports capacity-weighted mean, minimum and
maximum solid temperatures alongside actual fluid heat and combined balance.

## Spatial finite-wall conduction on native thermal slabs

`ThermalLbm::set_wall_conduction` and `ThermalLbm3::set_wall_conduction` prepare authored reciprocal conductances between existing finite wall cells. The original wall values remain the only physical solid temperature/storage owner. The shared sparse LDLT factor is prepared with diagonal scaling and explicit node/link/factor bounds; each subsequent boundary solve reuses retained workspace. Unconnected wall links retain their original scalar law.

A single backward-Euler boundary solve couples solid storage, sources, lateral or explicitly authored cross-plane conduction and the actual post-collision scalar reflection. Each edge is evaluated once with opposite endpoint contributions. Independently measured reflected heat, actual representable endpoint storage and internal conductive transfer must satisfy local and whole-plane first-law checks. Wall replacement invalidates the operator; failed preparation preserves the accepted plan. Checked source snapshots retain the plan and roll back wall state and transport history with the native populations.

This is fixed-capacity, fixed-resistance conduction on lumped wall cells. It does not resolve solid depth, infer geometry/conductivity, establish nonlinear material behavior, or remove the fluid solver's stability/refinement requirements. Time integration is first order for this solid solve.

## Finite solid cells on native thermal walls

`ThermalWall::FiniteSolid` couples an explicitly sized solid cell to each
normal D2Q5/D3Q7 slab link. Capacity, bulk temperature, through-wall/contact
resistance and external heat per original lattice step are authored in the
same scalar/energy coordinates as the fluid. A uniform declaration expands
once into the existing per-link spatial state; live bulk temperatures are
read through `thermal_wall_fields`. Ordinary steps allocate no wall scratch.

The local backward-Euler solve closes both `C dT = Q - q` and
`T_bulk - T_surface = R q` against the original reflected-pair heat. The
fluid receives that heat once. `ThermalWallExchange::solid` retains independently
measured solid storage, its external source, the actual fluid exchange and
the largest per-link first-law defect. Checked source/material steps admit
surface and solid-bulk temperatures and both summed and local balance before
publishing; existing spatial snapshots roll wall state back with flow, scalar
populations and advective history on any endpoint refusal. At extreme capacity
ratios, an unrepresentable solid temperature increment can leave a first-law
defect. Checked source/material steps refuse defects above the caller's balance
tolerance; ordinary unchecked steps expose the defect in their wall report.

Finite-wall kernel tests exercise the analytic interface/storage identities,
heating/cooling, conductivity changes and tiny/large capacity limits.
The shared scalar integration regression evolves both lattices under tangential
advection and spatial heating/cooling, checking fluid/solid heat accounts each
step and the independent capacity-weighted equilibrium after the source stops. Native
`finite_thermal_walls` regressions cover both lattice stencils, capacity-weighted
equilibrium, combined energy, one-time source deposition, independent wall
cells and endpoint rollback; the source-history test also rejects after full
candidate advancement. This is lumped first-order solid storage behind a
stationary planar resistance: no lateral wall conduction, resolved solid
diffusion, temperature-dependent solid law or coupled second-order time claim.



## Midpoint material tensor diffusion

The D2Q5/D3Q7 material slabs expose `step_with_material_tensor_diffusion_midpoint`
for source-bound, rotated K(T). The original native step is solved at the actual
old/new temperature midpoint with bounded under-relaxed iteration. Temperature
consistency and a fresh relative check of every principal conductivity must pass;
a strong axis cannot hide a mismatch on a weak one. The original caloric storage,
physical source/force conversion, wall exchange and energy acceptance remain
authoritative. No directional average replaces the tensor.

Each trial uses the original source-step transaction and accepted populations,
including force frames and advective history. Nonconvergence, source-domain or
endpoint refusal, and cancellation publish no partial time step. Density/Cp remain
constant over their source domains. This closes conductivity lag in the stated
fixed-point equation; it does not claim general second-order time accuracy or
unconditional stability of the split LBM. G1/G3/G4 integration cases compare an
independent nonlinear two-cell tensor reference, exact constant-law replay in
both dimensions, weak-axis admission, and transactional failure/replay.

## Planar fluid force-and-couple impulses

`Grid::apply_wrench_impulse` and `FreeSurface::apply_wrench_impulse` accept a total
linear impulse and out-of-plane angular impulse on strictly increasing original
wet-cell indices. Actual population mass, or retained partial liquid mass at an
interface, determines the centroid and scalar moment of inertia. The resulting
translation plus rotation is the least mass-weighted velocity increment realizing
the wrench. Cell centers are `(x+1/2,y+1/2)` in lattice lengths; native planar
impulse and work are per unit depth. A single cell cannot carry a free couple.
Periodic topology refuses because no unwrapped angular coordinate chart is owned.

The shared population impulse retains the original Guo force frame, intrinsic
density and physical viscous stress. Independent linear/angular momentum and
signed kinetic-work readback must meet separate budgets before publication;
invalid inertia, excessive speed, unsupported cells and cancellation preserve
accepted state. Partial liquid mass, fill and classification are retained exactly.
The implementation uses existing population scratch, one compact stencil and a
scalar inertia solve, without a dense impulse field or additional collision.
No time, heat or reciprocal body reaction is assigned. Six focused integration
tests cover the physical quadratures, partial inertia, stress, origin covariance,
Galilean work, braking, unresolved couples and transactional refusal/zero input.

## Compact borrowed free-surface impulse staging

`FreeSurface3::prepare_impulse_stencil` retains a mutable borrow of the original native owner while sampling a bounded, strictly increasing wet-cell stencil. Native populations, located retained liquid mass, physical Guo-frame velocity and original cell centers feed the shared population impulse kernel. Gas, inactive/static solid, piston, unresolved interface inventory, excessive speed and unrepresentable states refuse. Staging preserves original populations until the final commit; dropping a stage/candidate or final cancellation leaves the original unchanged. Candidate readback exposes each original mass/position and actual endpoint velocity so a finite receiving owner can independently gate reciprocal angular momentum and work before joint publication. No fluid clock, held force, fill/classification, geometry or unlocated carry is changed.

`MaterialFreeSurface3::prepare_native_impulse_stencil` reuses original physical impulse admission and scenario speed, while explicitly returning native-unit staging for a composing owner’s SI transaction. It does not silently reinterpret lattice readback as physical units.

The native stencil's three focused tests execute actual retained populations, force-frame stress, readback and transactional refusal. A further native refinement regression checks retained solid keys outside the physical domain and beyond representable child Morton coordinates: refinement refuses before allocation, preserving source mass, geometry and populations. All six native refinement cases pass after this guard. The prior three SI receipt checks remain unchanged; full streaming/material integration was not executed.

## Temperature-activated molecular wetting

`MolecularKineticLaw::try_new_arrhenius` optionally evolves the equilibrium jump
frequency as `K(T)=K_ref exp[(E_a/k_B)(1/T_ref-1/T)]` inside the original calibrated
temperature domain. Activation energy is explicit joules per molecular event;
the frequency in `MolecularKineticControls` is its value at the declared reference
temperature. This constant-prefactor Arrhenius closure adds no inferred adhesion
measurement, Eyring prefactor, extrapolation or activation-heat source. The
original constructor still selects constant frequency. Reference-temperature and
zero-activation responses preserve that original model exactly.

The same law feeds the existing dense, sparse and enthalpy wall-angle owners,
retaining original interface claims and actual local temperature. Accepted samples
report the local jump frequency separately from static material measurements.
The rate and velocity scales must be finite and positive across both temperature
bounds; thermal domain refusal retains native flow, enthalpy and wetting history.
Per-length molecular loss remains diagnostic, without claiming a resolved line
measure, a global surface-energy balance or heat deposition. G1/G3 unit cases
check independent forward kinetics, low-speed friction and constant-model limits;
the existing thermal integration case also exercises activated rates and rollback.

## Isothermal kinetics from reversible equilibrium adhesion

The dense D2Q9 and sparse D3Q19 material owners can attach an original equilibrium
angle or Young--Dupre work with `with_equilibrium_planar_wetting`, then explicitly
select Cox or molecular kinetics with `with_dynamic_equilibrium_wetting`. The
same frozen liquid tension feeds both adhesion and native pressure. Static
coefficients sampled at fixed external conditions retain their source axes and
receipts; time, loading and rate restrictions cannot be silently frozen. No
advancing/receding measurements or contact-speed query are fabricated from work.

`with_equilibrium_wetting_hysteresis` optionally binds independent static pinning
thresholds from that interface, preserving the equilibrium potential and separate
resistance receipt. The existing full-angle rate domains, geometric speed history
and transactional native publication still apply. Pure equilibrium alone bypasses
the velocity quotient when fill changes on a locally flat wall stencil; adding
kinetics or pinning retains the unresolved-speed refusal. The backward geometric
speed drives the next step's capillary ghosts. Line losses remain diagnostics;
these owners do not infer a line measure, heat deposition or closed energy balance.

Four public integration groups exercise both source origins and rate laws in
both dimensions, independent forward equations and native response, separately
sourced pinning, domain rollback, and the pressure-loaded flat-fill distinction.

## Heated wetting at neutral transverse corners

The height-graph equilibrium-wetting owner now closes the two ends of its
sourced side wall against the explicitly selected neutral Cartesian transverse
walls. The simultaneous conditions `n.n_wall=cos(theta)` and
`n.n_transverse=0` give zero transverse height slope and the existing sourced
normal slope. Graph normals, area, Marangoni projection and public heat patches
share this metric. The original raw/smoothed VOF pressure ghost uses even
reflection at that same neutral wall, including diagonal corner samples.
The internal boundary selection survives thermal angle-field replacement;
ordinary planar wetting retains its original stencil.

This closure requires one sourced stationary wall, no internal solids or piston,
complete monotone liquid/gas columns, one angle through each occupied interface
band, and the original graph-slope envelope. It does not supply an angle for a
second substrate, arbitrary wedges, corner films, corner pinning, capillary heat
or a closed surface-energy balance. Invalid source work, geometry and cancellation
still refuse the complete thermal/momentum candidate. G1 tests cover both
graph orientations on all Cartesian wall pairs, raw affine-plane preservation,
zero artificial Laplace pressure, independent corner traction and an actual
heated Young--Dupre corner response through the public material owner.

## Material-flow enthalpy reconstruction

`MaterialFlowAdvection::with_reconstruction` attaches the original thermal owner's shared `LimitedEnthalpyReconstruction` to solved material face fluxes. Selected contiguous axial cells reconstruct total specific enthalpy, including a declared passive equal-density latent inventory; paired face transfers, source-law inversion, accepted history and the combined conduction/advection timestep restriction remain with `fs-thermal`. Unselected cells and external streams retain donor enthalpy. Empty plans retain the donor-cell arithmetic. The observed timeline retains first-order time integration on stationary control volumes; it does not infer multidimensional junction slopes, phase mobility, density-change flow, pressure work or moving interfaces. The `advection()` accessor exposes the underlying graph and actual transport audit; implicit clients must attach the same `reconstruction()` plan to their own native implicit transport owner.

Pure frozen `MaterialFlowAdvection<_, ()>` owners additionally expose native
`advance_ssprk2_to` and `advance_ssprk2_adaptive_to`, including corresponding
`with_coefficients` entries for repeatable stage sources and link/bath
conductances. These methods retain the solved `rho * face_volume_flux` graph,
original material/density snapshot and selected reconstruction in the returned
accepted wrapper. Original fluid-cell mass, geometry and caloric-law identity
are admitted at stage queries and before publication. The native thermal owner
alone integrates storage and heat/mass accounts. Adaptive comparisons are one
full step versus two fine steps, and publish only the conservative fine state
and its ledgers. Mandatory datum-invariant h error resolves latent storage at
unchanged temperature; optional T error queries the original caloric law.
Native events, CFL/caloric admission, shared work budgets and cancellation
remain authoritative. The smooth-solution local estimate gives no spatial or
certified global accuracy bound. Stateful inner observers and live momentum
feedback are not repeatable coefficient queries and do not expose these methods.
Integration regressions use actual material-driven momentum for analytic
latent-throughflow temporal refinement, reconstruction in both flow directions,
separate supplied/stream heat, map/law admission and transactional refusal.

`MaterialFlowFeedback::with_reconstruction` retains that same shared plan across
live source viscosity, bounded Boussinesq buoyancy, mushy resistance and retained
momentum updates, including both native implicit entries. Each actual source
epoch rebuilds the original mass graph before binding the same reconstruction;
no additional momentum solve or alternate caloric law is introduced. Required
axial interior connections remain present at zero flow, allowing rest and flow
reversal. They count against the original connection budget. Dormant transverse
faces add no artificial junction, while nonzero crossflow remains visible and
refuses an unsupported axial stencil. A missing neighbor at an exact-solid lock
also refuses the candidate; external or cross-solid connectivity is not inferred.

The stationary-skeleton liquid-carrier attachment reconstructs phase-carried
liquid enthalpy, leaving solid latent inventory in original mixture storage.
Its equal-density, pure-substance, source-support and physical Stokes-work gates
still apply. The optional bulk entropy face limiter remains incompatible with
liquid-only carriage. Reconstruction can be selected before the first advance;
accepted feedback preserves it with the original material pins, physical history
and transactional candidate. Empty selection preserves donor arithmetic and CFL.

## Adaptive implicit thermoviscous flow

`MaterialFlowFeedback::advance_implicit_adaptive_to` controls time error in the
existing quasisteady Stokes/material/enthalpy endpoint solve. Each comparison
uses one original backward-Euler interval and two fine intervals; only the fine
thermal states and their actual endpoint flow, source and physical work are
accepted. The original `ViscousHeating` owner integrates both accepted intervals,
retaining wall/pressure/body supply and molecular/Darcy/Forchheimer heat without
turning momentum residuals or temporal error into heat. All coarse, fine,
rejected and finite-difference evaluations share the blueprint's one
`max_solves` budget. Individual one-step reports retain per-interval work counts.

Mandatory enthalpy error scales changes from the comparison's starting state,
so latent inventory remains resolved at fixed temperature and absolute
formation-enthalpy datums do not inflate the allowance. Optional temperature
error queries the original caloric laws. First-order step-doubling uses the
actual representable fine spans and accepts no extrapolated state. The caller's
whole-advance joule allowance bounds cumulative absolute accepted algebraic
residuals and the final compensated, datum-relative storage/heat balance;
native actuator/physical-loss power tolerance remains separate. Reports retain
the final original coupled endpoint and accumulated accepted physical accounts.

This entry retains fixed topology, stateless `()` material feedback, and the
original stationary-heating restrictions. Interior or unprocessed native events
refuse before momentum work; target-time events remain pending. The native
accepted-substeps budget counts both fine intervals. Minimum timestep, trial
budget, source/work/nonlinear refusal and cancellation publish no partial
thermal or flow prefix. The estimate controls local time error on admitted
smooth solutions, not spatial/model error or a certified global bound. Focused
integration cases compare continuous source-viscosity heating, latent mushy
mobility under shifted h datums, and shared speculative-work/transaction gates.

## Hydrostatic rest under material density-dependent gravity

`MaterialFreeSurfaceScaling::free_surface3_hydrostatic` initializes a planar
column using the original material-bound Guo force, collision, streaming and
gas-pressure reconstruction. Gravity points toward one declared Cartesian
supporting wall; zero gravity is also admitted. Complete liquid layers retain
the original stationary solid tiles and at least one layer of gas. The caller
supplies a positive absolute pressure reference and a density-excursion bound.
Capillarity and wetting must be explicitly absent.

`free_surface2_hydrostatic` applies the same original density recurrence to
the dense D2Q9 owner. Its X/Y supporting wall is halfway between the wall row
and first interior row; the declared depth counts only complete interior
liquid layers, with at least one interior gas layer. Boundary interface cells
retain full mass so their volume agrees with the native gas-pressure face.
This explicit constructor does not change the generic classifier's half-filled
interface convention. Observations and mass remain per unit out-of-plane depth.

`FreeSurface::new_hydrostatic_column` applies that same profile to a fresh
native D2Q9 grid with an existing wall/weir mask. Complete occupied columns
run from the floor to one declared upper gas face; vertical/submerged gas
exposure refuses. The input is closed and nonperiodic, with downward held
gravity, zero external force and raw unit-density BGK rest populations.
The constructor retains the mask, initializes the original physical rest
frame and full interface mass, and introduces no per-step work. Its masked
weir regressions share the dense integration seam below. Historical generic
initialization and its half-filled inventory remain unchanged.

For lattice acceleration magnitude `a`, set `h = a/(2*cs^2)`. The nearest wet
layer has `rho = 1/(1-h)`; successive deeper layers multiply by
`(1+h)/(1-h)`. This surface factor follows the native missing-gas-link rule.
A continuum pressure profile evaluated at cell centers alone does not supply
the same discrete rest state. Layer densities are computed once and reused
while the original constructor owns liquid mass, sparse activation and the
physical-zero-velocity Guo frame. Ordinary `step()` retains the same gas datum.

This balances the weakly compressible numerical operator. It is not a liquid
equation of state or an incompressible mass projection: native mass contains
the retained density perturbation. Pressure and submerged traction approach
their incompressible hydrostatic limits as `dt/dx` decreases with the physical
configuration fixed; refining `dx` at fixed `dt/dx` does not remove that
compressibility error. The constructor adds no moving-solid or partially
wetted-body force closure.

G1/G3/G4 regressions in `tests/material_free_surface3/hydrostatic.rs` exercise
all six gravity directions, zero-gravity equality, pressure-datum invariance,
actual reflected-link traction and torque with decreasing compressibility,
mass retention and transactional refusal. The declared density, source-domain,
geometry and cancellation bounds remain with the original owners.

The dense `tests/material_free_surface/gas_pressure2.rs` cases exercise all
four directions including a single wet layer, full-depth mass retention,
ordinary native rest evolution, the `dt/dx` pressure convergence ladder,
zero gravity, density/viscosity similarity, pressure-datum invariance and
invalid/cancelled construction. Both dimensions use one temporary profile
value per wet layer and add no work or retained storage to time stepping.

## Conservative streamwise laminar-film evolution

With `thin-film-evolution`, `film::evolution::ThinFilmLine` advances actual
strictly wet cell thickness by `h_t + q_x = 0` on a uniform
line. The face force combines held tangential gravity with the gradient of
`p = rho*g_normal*h - sigma*h_xx`. Original Newtonian `ShearedFilm` or prepared
thermal Carreau-Yasuda slit profiles supply the flux, surface speed, wall
stress and physical viscous heat together. Each original Carreau face source
is prepared once at its frozen temperature; changing height reuses that
preparation. Pressure forcing and surface traction enter one nonlinear
profile without apparent-viscosity replacement or superposed flows.

Godunov selection evaluates both neighboring heights and, when present, the
original flux's interior turning point `h = -tau_surface/force`. This admits
opposing body/shear forcing and counterflow without assuming a monotone flux.
The same selected profile supplies all work accounts. Matrix-free backward
Euler solves normalized height changes with the native fallible finite-
difference/Newton-Krylov owner. Every face transfers equal and opposite
volume between cells; closed ends have no flux and natural zero-slope
capillary boundaries. No clipping, mass repair or post-solve rescaling enters.

The step independently admits each cell balance, total volume change and
`Delta E + E(delta h) + viscous_heat - body_work - surface_work`, all in
physical SI units per unit out-of-plane width. `E` is the discrete quadratic
gravitational plus capillary energy. The nonnegative `E(delta h)` is the
backward-Euler numerical energy defect and remains separate from actual
viscous heat. These accounts use the accepted native profiles, not rejected
Newton or Godunov candidates.

Authored wetness, height/physical-length and slope limits screen the accepted
opening state and every proposed height. Candidate profiles also pass a
reduced-inertia limit using actual profile extrema, including interior
counterflow, and the minimum original viscosity; its ceiling must be below one.
Height, slope and
inertia screens are model admission, not calibrated accuracy guarantees.
Whole-call endpoint/profile budgets include all finite-difference probes and
discarded face states; native constitutive work remains separately bounded.
Cancellation, unsupported sources, failed solves or failed conservation gates
preserve all accepted heights and time, including the final publication poll.

`MaterialFilmLiquid::thin_film_line` retains the exact borrowed original liquid
source and its receipts while passing its original SI density, Newtonian
viscosity and surface tension into this same owner. The reduction assumes a
stationary planar wall, quasi-steady wall-normal momentum, small slopes,
fixed density and uniform surface tension. It omits startup momentum storage,
wetting/contact lines, rupture, evaporation, transverse structure and thermal
evolution. Face temperatures and imposed surface traction are authored held
inputs. First-order time/transport still require refinement for accuracy.

Five focused production-source regression groups pass: independent discrete
Fourier gravity/capillary decay, closed leveling and energy, original shifted
Carreau profiles against independent quadrature, both signs of the interior
Godunov turning point, and atomic source/budget/final-cancellation retry.
The numerical source-linked build executes the original material, math and
solver implementations; runtime-context paths are uninhabitable in that
projection. The source-bound liquid integration regression is added but its
full native crate/dependency execution is not claimed.

### Open film connections with original viscous losses

`ThinFilmBoundary::Open` and `ThinFilmLine::advance_open` admit two explicitly
held hydraulic pressure potentials and surface tractions. Each physical end
uses the adjacent accepted thickness and the original profile over half a cell:
`F_left = rho*g_t - (p_first-P_left)/(dx/2)` and
`F_right = rho*g_t - (P_right-p_last)/(dx/2)`. Natural capillary end slopes
remain zero. Boundary pressures use the same gravity/capillary potential datum
as the cell energy; they are not arbitrary gauge offsets or inferred inputs.

These two half-cell connections complete the physical domain length and carry
their own original viscous heat, body work and surface work. Their signed flows
enter the boundary cells' finite-volume equations. The independent volume gate
is `Delta V - dt*(q_left-q_right)`. The energy gate additionally subtracts the
actual reservoir work `dt*(P_left*q_left-P_right*q_right)`. Positive flow is +x
at both ends, so either end may supply or remove liquid. No prescribed-flux
shortcut, unmodeled boundary resistance or residual-to-heat conversion enters.

`ThinFilmProfile` is the common physical response for interior and boundary
connections. `ThinFilmFace` retains its from/to cells and this profile;
`ThinFilmOpenReport` retains both original end profiles and their authored
pressures/tractions. Volume change and volume balance residual remain distinct.
Each original Carreau source is still prepared once: an open N-cell line needs
N+1 authored temperatures ordered left end, N-1 interior faces, right end.
End temperatures are never inferred from an interior source. `face_count()`
counts interior links and `profile_source_count()` includes the two open ends.

The same source support, profile-extremum inertia admission, bounded numerical
work and final cancellation gate apply to both ends. Missing/invalid end data,
unsupported sources, unresolved equations or failed accounts preserve accepted
heights/time. `MaterialThinFilmLine::advance_open` delegates through the same
original liquid binding and preserves its source receipts. The open extension
retains first-order quasi-steady long-wave scope; it does not resolve wetting,
contact lines, entrance development or the downstream free jet.

Focused open-film regressions cover uniform forward/reverse throughflow over
the complete domain length, independently manufactured nonuniform endpoints
with pressure work and capillarity, original nonlinear Carreau end profiles,
and atomic source/control/cancellation refusal. They accompany the existing
closed/periodic physical cases without changing their acceptance tolerances.
All nine evolution cases passed in the existing build of original numerical
sources; no compiler warnings were emitted. This does not establish the full
asynchronous runtime or original liquid-source adapter integration.

## Original-law shear-free Carreau film base

`film::carreau::CarreauFilm` resolves a positive prescribed flux per width at
fixed temperature under tangential gravity. It owns the original thermal
Carreau-Yasuda law, source-rate ceiling and accepted `CarreauSlitProfile`.
The lower half of a symmetric stationary-wall slit of gap `2h` supplies the
shear-free film exactly: `tau(y)=rho*g_t*(h-y)`. Native constitutive inversion,
positive quadrature, temperature shifting, source support and viscous heat
remain authoritative. No fitted effective viscosity enters the profile.

A bounded thickness solve brackets against the original zero-shear Nusselt
thickness and supported wall stress. Safeguarded Newton uses the exact
`dq/dh=h*gamma_wall` sensitivity and admits both flux residual and relative
thickness correction. Requested and accepted flux remain distinct. Point
samples retain apparent and differential viscosities at the same local shear;
`mu_t*u''=-rho*g_t` supplies curvature. Surface shear is zero, wall speed is
zero, and original heat per wall area equals gravitational work within the
native quadrature tolerance. Density and surface tension are explicit actual
SI inputs, not inferred from a constitutive curve or fabricated material card.

This fully developed, isothermal, planar base does not itself claim stability,
entry-region accuracy, thermal feedback or a generalized turbulent threshold.
All support/work failures and cancellation refuse without a partial film.
Focused regressions compare the exact constant-source Nusselt limit at shifted
temperature, an independently integrated nonlinear profile, gravity/heat balance,
source limits and cancellation.

## Lagrangian turbulence with source-resolved rheology

`RheologySgs3::Lagrangian` composes the original source-resolved molecular
viscosity with pathline-averaged dynamic SGS on periodic or stationary channel
grids. `ChannelRheologySgs3::Lagrangian` uses the same history with prescribed
moving channel walls. Source shear, turbulence fitting and pathline departures
use the same incoming Guo force frame and physical wall velocities. Material
queries retain SI shear and time; the authored SGS memory time remains in
lattice timesteps. Each step samples molecular transport afresh, without
feeding an earlier augmented viscosity back into the source law.

The original source and filter workspaces own all numerical storage. The
history trial remains unpublished through collision, streaming, endpoint-force
admission and moving-wall work measurement. Only a successful native step
commits its populations, force/clock, source report and turbulence memory;
rejected candidates preserve the elapsed history interval for a retry. The
`Lagrangian` report variant distinguishes this model from instantaneous fits.

Focused G0/G3/G4 regressions compare several evolving periodic/channel steps
with the independently composed native operators, serial/pooled results,
source/force refusal and retry, and thermal Carreau--Yasuda moving channels
in all three orientations. This retains the exploratory closure's topology,
coefficient, molecular-ratio and departure bounds. It does not supply
unresolved non-Newtonian stress correlations, a thermal energy equation or
turbulent-DNS validation. Locally adaptive memory is an explicit option below.

## Physical strain in planar material rheology

Native and source-backed D2Q9 thermal and free-surface rheology use the same physical
`sqrt(2 dev_3(D):dev_3(D))` invariant as `fs-material`, embedding the original
planar tensor with `D_zz = D_xz = D_yz = 0`. Equal in-plane expansion therefore
retains its actual three-dimensional deviator. Guo half-source correction and
independent central-moment trace relaxation precede this embedding; the latter's
planar trace contributes `div(u)^2/3` to the squared material shear rate.

Thermal power separately retains `2 mu dev_3(D):dev_3(D) + zeta div(u)^2`.
Ordinary tied collision keeps its existing implicit `zeta = 2 mu / 3`, so its
total planar loss remains `2 mu D:D`; the material query no longer counts that
bulk term as shear. Independent planar bulk relaxation still requires a physical
bulk-heating owner and is refused by the shear-only heating entry. Native
lattice-unit laws and SI-converted `PhysicalRheology` use this same material
invariant, so a planar flow does not acquire another apparent viscosity when
embedded in D3Q19. This corrects the former native full-strain convention;
incompressible simple shear and Poiseuille retain their original law.
Manufactured original-population tests cover native and source
owners, force correction, trace relaxation, actual viscosity, and heat partition
against the material tensor and independently prescribed thermal source.

## Density-aware coarse/fine interpolation

`RefinedChannel` interpolates spatial and temporal ghost data in the physical
Guo population frame `h=f+S(u,F)/2`, whose first moment is `rho*u`.
Each donor's actual density and force determine its source transform; the
receiving force frame is applied only after the physical interpolation.
Restriction likewise blends rescaled physical populations before rebinding
the coarse force. This preserves the existing mass-weighted velocity and
physical second moments without a force/velocity covariance artifact at
nonuniform density. Unforced transfer arithmetic is unchanged. Existing
rescaling, reflux, timestep order and low-Mach model limits still apply;
this does not introduce a variable-density incompressible equation of state.

## Spatial pressure reservoirs and sealed D2Q9 channels

`core2::PressureProfileX2` borrows one positive finite lattice-density value
per y row at each x face. `Grid::step_pressure_profile_x` and its selected-wall
receipt variant apply these held pressure profiles through the same native
regularized reconstruction as `PressureX2`. No callback runs inside a step and
there is no profile allocation. Uniform profiles reproduce the uniform-pressure
arithmetic. All candidates reuse collision scratch; its O(ny) extra capacity
is retained while the caller-visible extent returns to the collision field.

Pressure faces now admit paired stationary-wall face/interior rows. Non-periodic
y requires complete top and bottom wall rows, and each x face needs at least
one open fluid row. Velocity/pressure inlet mode retains periodic y and fluid
face requirements. A profile-length, topology, nonfinite force, or late donor
refusal leaves grid populations unchanged. Wall populations are never rebuilt
as pressure cells. Wall receipts retain the original bounce-back link impulse;
the subsequent reservoir reconstruction is a separate momentum exchange.

Each fluid face prescribes `p = cs² rho`, copies its streamed first-interior
physical velocity and non-equilibrium stress, and restores the destination's
own `rho*g + fext` half-source. The pressure profile can therefore sustain
transverse hydrostatics with the discrete relation
`rho[y+1]/rho[y] = (1 + g_y/(2 cs²))/(1 - g_y/(2 cs²))`.
This follows directly by equating opposite Guo half-source populations across
a streaming link; the continuum exponential alone is not the discrete balance.
Pressure remains a weakly compressible isothermal lattice quantity. SI loads
require a shared material pressure reference. These regularized reservoirs
reflect acoustics and do not supply a characteristic or thermodynamic boundary.

The formulation retains the equilibrium/non-equilibrium split of Guo, Zheng
and Shi, [Chinese Physics 11 (2002), 366–374](https://doi.org/10.1088/1009-1963/11/4/310),
with the existing second-order Hermite stress projection of Latt and Chopard,
[Mathematics and Computers in Simulation 72 (2006), 165–168](https://doi.org/10.1016/j.matcom.2006.05.017).
Focused `pressure_driven2` regressions cover both signs of hydrostatic gravity,
the plane Poiseuille pressure/viscosity law, local Guo/stress recovery,
uniform-profile replay, native wall receipts, scratch reuse and atomic refusal.

## Original-law nonlinear tensor SDIRK2

`advance_nonlinear_tensor_conduction_sdirk2` advances the existing full-tensor
voxel conduction equations with two physical SDIRK2 stages. Every nonlinear
trial reevaluates the original conductivity interval integrals and endpoint
Newton tensors. Both stages retain the physical starting temperature in
storage; stage two uses the accepted first-stage temperature only as its seed.
Full anisotropy, finite-contact traces, sinks and independently authored
stage source/boundary samples retain their original spatial treatment.

The shared interval ledger integrates actual external heat with the tableau
weights. Internal stage heat is excluded from generated power, and prescribed
reservoir storage is counted once, including changing prescribed masks.
Krylov iterations, Newton corrections and rejected line-search candidates
each have one budget shared across the complete interval. Refusal or
cancellation returns no partial interval and leaves inputs unchanged.

Each cell has a fixed, complete sensible capacity C in J/K, without a second
volume factor. This entry does not supply temperature-dependent caloric
storage or latent heat. Second-order time accuracy requires smooth laws and
correctly sampled smooth forcing with fixed constraints; forcing jumps must
bound intervals. The L-stable tableau does not establish unconditional
nonlinear convergence, arbitrary-step positivity or a discrete maximum
principle for general rotated anisotropy.

## Locally adaptive Lagrangian momentum memory

`sgs3::dynamic::lagrangian::AdaptiveLagrangianMemory` is an explicit opt-in
configuration for the existing signed Germano history. Before its first
accepted fit, `LagrangianDynamicSgs::with_adaptive_memory` and the native
`LagrangianSgsWorkspace3` can author a positive dimensionless coefficient
and a minimum time. The step configuration's original `memory_time` is the
upper bound. Existing owners without this selection retain their fixed-time
arithmetic. Scalar and component owners expose the same momentum selection;
their independently authored heat and species memory times remain unchanged.

The local law is `T = theta Delta (I_LM I_MM / rho^4)^(-1/8)`, based on
Meneveau, Lund and Cabot (1996), DOI `10.1017/S0022112096007379`.
It uses the advected accepted history and the original grid-filter width.
The density normalization removes the squared-density units of these Favre
contractions; this variable-density extension is not a calibrated claim.
Both time bounds and theta remain explicit. Zero, undetermined or nonpositive
signed numerator history selects the upper bound without discarding its sign.
This bounded policy differs from the paper's positivity mechanism.

The first fit remains the original local spatial average. Each later local
time is frozen during the existing first-order pathline/source relaxation.
Logarithmic invariant arithmetic avoids products, fourth powers and premature
underflow; old-history decay enters the shared scale before exponentiation.
No per-cell arrays are added. Accepted minimum/maximum times publish only
when the original turbulence trial commits, including native collision and
material-rheology consumers. Cancelled, dropped or refused trials preserve
accepted histories and times. Native-grid times remain in lattice timesteps.

Six new regressions verify unit similarity, advected-history relaxation,
time bounds, signed and extreme-scale histories, all channel orientations,
transactional retries and independent thermal/species memories. The exact
production SGS module graph passes all 80 tests with the pinned compiler and
actual native math dependencies. The full fs-lbm grid/material adapters were
source-reviewed but not compiled in the constrained environment. This adds no
adaptive scalar memory, backscatter-energy equation, calibrated theta, DNS
qualification or higher-order pathline integration.

## Source-local planar bulk viscosity

The D2Q9 thermal material owner now admits separately authored physical
`mu(T,gamma)` and `zeta_3(T,gamma)` through the existing `BulkViscousMaterial`.
Planar reduction uses `zeta_2=zeta_3+mu/3`, so even constant physical bulk
viscosity requires a local trace relaxation when shear viscosity varies.
The original central collision reads an optional prepared per-cell trace-tau
field; default/native scalar selections retain their original arithmetic.

The bounded material solve uses the original Guo-corrected stress numerators
`A=gamma_dev2*tau_s`, `B=sqrt(2/3)|div u|*tau_trace`, and
`N=hypot(A,B/sqrt(2))`. Its auxiliary relaxation satisfies
`gamma=N/(lambda*dt)` and
`1/lambda=hypot((A/N)/tau_s,(B/(sqrt(2)*N))/tau_trace)`.
Normalized weights avoid squared-component overflow/underflow; zero stress
queries the zero-rate source. The report's relaxation defect is in lambda.
Collision receives the actual shear and trace taus from that same source,
never lambda. Shear-only and trace-only inputs recover their respective taus.

The unchanged thermal transaction converts `mu*gamma²+zeta_3*(div u)²` once,
retains original material receipts, checks the endpoint, and publishes the
applied trace field with the accepted state. Explicit `zeta_3=0` is allowed
in planar flow through the positive `mu/3` trace coefficient. The 3-D trace
collision retains its strict positive-relaxation admission; negative bulk
viscosity and unrepresentable positive losses refuse. This remains explicit,
fixed-density/Cp low-Mach sensible-heat transport, without a total-energy claim.
Manufactured tests cover the physical tensor/heat identity, source-dependent
coefficients, zero bulk, original native collision, and full rollback.

## Lagrangian turbulent heat with original material sources

`thermal::sources::sgs3::lagrangian` composes the shared joint Favre/Germano
momentum and scalar histories with the original D3Q19/D3Q7 source transaction.
The resolved slab retains periodic tangent directions and its actual rigid
halfway walls. Current prescribed-temperature, insulating, flux, film and
spatial wall derivative policies are shared with the instantaneous channel
closure; the native thermal owner alone imposes and measures their wall heat.

Stress and heat keep independent signed contractions, memory times, coefficient
limits and molecular-ratio limits. A zero momentum fit cannot disable a
positive heat fit. The existing optional invariant momentum memory is exposed
without altering scalar memory. Sampling reads one accepted Guo-centered
density/velocity/temperature state before replacement forces. The caller
supplies actual sampling time in lattice timesteps; ThermalLbm3 has no separate
authoritative clock. Elapsed sample time controls bounded pathline departures.

One prepared token borrows only the reusable workspace and retains both
candidate histories plus admitted total relaxation fields. Native source, heat
balance, temperature range, endpoint Mach and cancellation checks precede
publication. Material wrappers additionally retain that token through every
original source epoch and outer publication check. Only successful completion
is followed by infallible joint history commit. Dropped, cancelled and refused
steps preserve both histories, including their last sample time for retry.

Material APIs cover ordinary source-resolved transport, source-density
Boussinesq forcing and endpoint-averaged source-density forcing. The last keeps
all three original material epochs: initial transport, thermal-predictor force
quadrature and the actual accepted endpoint. Original SI source samples, fixed
reference inertia/sensible storage, density/Cp bounds and force conversions
remain intact. Fresh molecular mu/rho0 and k/(rho0 Cp0) own their ratio limits;
augmented rates never feed back into a later molecular source query.

Field arrays are prepared once under an explicit numerical byte budget and
reused by sampling and history publication. Original native snapshots retain
their independent cell budget. The instantaneous and Lagrangian adapters share
the same rate augmentation and shear-power arithmetic. No alternate collision,
wall solver or material law is introduced. This remains explicit, low-Mach
sensible-temperature transport; modeled kinetic transfer is diagnostic rather
than heat deposition, and no total-energy or turbulent-DNS claim is added.

Five new native regressions compare multiple evolving steps against the
independent shared history owner and original transport, including fixed and
adaptive momentum memory, independently positive heat, actual wall policies,
Guo frames, exact zero controls and late refusal/replay. Three material groups
cover all wrappers, fresh source epochs and endpoint buoyancy, then warmed
history rollback after late source refusal and cancellation. All changed Rust
files passed pinned parser checks and bounded source review. These eight new
integration regressions have not been compiled or executed because required
dependency artifacts and build storage are unavailable; earlier shared-kernel
and atmospheric test results do not establish this native integration gate.

## Physical planar liquid wrench impulses

`MaterialFreeSurface2::apply_wrench_impulse` applies an explicitly authored
linear impulse and independent out-of-plane couple to retained liquid cells.
Inputs and receipts use the existing per-unit-depth convention: linear
impulse in kg/s, angular impulse in kg m/s, liquid mass in kg/m, and signed
work in J/m. Positive angular impulse is counterclockwise. No thickness or
extra timestep is inferred. The original source density, cell length and
actual interface fill determine the selected mass, centroid and inertia;
the native minimum-mass-weighted translation plus rotation realizes the wrench.

Conversion transports the requested moment to the representable native
origin and transports actual deposited and observed moments back to the
original authored SI origin. Separate physical linear, angular and work
budgets admit the complete SI receipt before native population publication.
Scenario speed, topology, singular-couple and cancellation gates retain the
original transaction; a late unit or balance refusal also leaves the accepted
liquid unchanged. Held Guo forcing, viscous stress, density, fill and material
identity remain with their original owners. No collision, thermal advance or
reciprocal body reaction is implicit.

Focused regressions cover independent physical impulse/inertia/work sums for
bulk and half-filled cells across distinct length/time scales, covariance
under a shifted SI origin, signed braking work and state-preserving native
and physical refusals. Dependent-receipt rejection is checked at the native
publication boundary.

## Original caloric storage in nonlinear tensor conduction

With `enthalpy-coupling`, `advance_caloric_tensor_conduction` and `advance_caloric_tensor_conduction_sdirk2` add backward Euler and SDIRK2 storage in the original caloric law's specific-energy coordinate. Each cell retains its supplied reference mass and original energy basis/datum. Newton solves centered energy increments with actual storage `m * (e - e0)`; face traces remain Kelvin. The original analytic `dT/de` scales only the cell column of the full-tensor Kirchhoff tangent. Exact zero slopes on latent plateaus remain regular through the positive mass/storage diagonal, without apparent-capacity smoothing, numerical differentiation, or a prescribed spatial phase front.

`TensorCaloricLaw` extends the canonical `ThermalLaw` with original bounded energy support, an analytic temperature derivative, a fused constitutive query, and reference-cell admission. Direct implementations retain the equilibrium phase curve, restricted pure-solid view, and original heat-storage assembly; amount-bound storage enforces its original mass normalization. Every prescribed-temperature cell requires an explicit original specific energy at each stage, validated against its authored temperature, so an isothermal phase fraction is never selected by a temperature inverse. Prescribed boundary-face temperatures do not require a cell-energy prescription.

Both stages reevaluate original caloric and tensor laws, share nonlinear/Krylov/backtracking budgets, and reuse the independently audited interval heat ledger with actual representable stored-energy changes and prescribed-cell reservoir corrections. Conductivity remains the independently supplied tensor law `K(T)`; independent phase-fraction dependence, moving geometry, expansion work, and advection are outside these entries. Backward Euler is first order. SDIRK2 is second order on smooth caloric/forcing segments with fixed constraints; unresolved sharp phase transitions or forcing changes require time refinement or interval splitting. Either stage's refusal or cancellation returns no partial interval.

## Lagrangian momentum and heat in conservative NSF3

`NsfGrid3::with_lagrangian_thermal_sgs` binds the original joint Favre--Germano
pathline histories to native D3Q27 mass/momentum and conservative total energy.
Periodic fluid boxes and resolved channels with two periodic axes are admitted;
actual translating-wall velocities and original thermal-wall derivative policies
feed the same filters. Independent momentum and temperature memories, optional
bounded invariant momentum times, coefficient bounds and molecular-ratio
admissions are explicit. Times are lattice timesteps from the native grid clock.
Source-owned viscosity, heat capacity, conductivity and bulk relaxation retain
their original state and domain. Only deviatoric collision receives modeled
viscosity; shared interior heat faces receive rho*(cp/R)*alpha_t. Mechanical
transfer stays in the original total-energy stress work, without an additional
heating source. Molecular half-cell conduction remains the wall heat law.

Every endpoint corrector/probe starts from the same accepted signed histories.
The final source/force/temperature state retains a borrowing joint history token
through all later native work, energy and cancellation admissions; one infallible
commit accompanies physical publication. Failed/retried steps cannot age memory.
Same-time wall/force changes re-emit the accepted contractions with zero elapsed
transport, without publishing a new model sample. Regular timesteps reuse the
264 nominal bytes/cell filter/history payload and do not copy history or repeat
the final fit. Thermal-wall preview intervals allocate an independent workspace
within their cell and original per-workspace bounds, copying only accepted
histories/scales/policy. A rejected multi-step prefix cannot modify its parent.

G1/G3/G4 regressions cover the one-update endpoint against the original shared
model, conservative closed-box energy, final-checkpoint rollback/replay,
moving-channel same-time refresh, independent multi-step thermal trials,
workspace/geometry/usage and representable-clock refusals. Core regressions
cover zero-time refresh with both fixed and adaptive momentum memory, all
three channel normals, current-strain re-emission and history-copy independence
from rejected scratch. An exactly zero interval is an identity; a positive
interval whose adaptive relaxation underflows still refuses.
The pathline/source update remains first order even though native energy uses
its original endpoint trapezoid. No open boundary, interior-solid filter,
backscatter-energy model, complete compressible LES closure, turbulent-wall
function, shock qualification or DNS claim is added. Native integration tests
are authored but require execution in a complete build environment.

## Original source contacts in nonlinear tensor conduction

The native nonlinear tensor owner consumes the existing ordered
SourceFaceContact assignments directly. This applies to steady conduction,
sensible backward Euler and SDIRK2, and both original-caloric time drivers.
Every finite source contact retains two physical face temperatures. The
original declared interface query seeds these traces; unequal bulk cell
temperatures may lie outside the interface source's temperature box, while
their conductivity intervals must remain supported by their original laws.

Every physical Newton trial refreshes the same selected original resistance
claim at those actual A/B surface temperatures. Card identity, surface order,
pressure, separation, finish, medium and other supplied conditions remain
with the existing source owner. Supported tangent payloads are a declared
constant scalar, an original linear curve, or original affine surface cells.
The analytical chain rule for q = A*(Ta-Tb)/R includes both resistance
temperature derivatives. Source corners use original one-sided primal
tangents and do not establish a unique adjoint derivative.

The two trace equations receive exact opposite contact heat and derivative
rows. Cell increments still condense locally into the existing sparse face
system; the contact adds no storage, external source, boundary heat or outer
fixed-resistance thermal solve. Numerical source responses and derivatives
are cached per trial, and source derivative preparation occurs once per
spatial solve. Original source receipts and actual bulk face heat are returned
in TensorConductionSolution.source_contact_states, with independent local
interface, contact and global energy checks.

TensorNonlinearConfig owns Newton, Krylov and rejected-trial work budgets,
shared across SDIRK2 stages. Scalar SourceContactControls govern only the
separate scalar outer iteration. Duplicate or overlapping contact owners,
unsupported original source states, failed convergence and cancellation
return no partial endpoint or interval. Frozen tensor entries continue to
require fixed contacts; existing fixed-resistance and perfect-interface
arithmetic retain their original paths.

## Original directional sources in tensor conduction

`RotatedSourceConductivity` borrows three original `ResolvedConductionLayer` owners in an explicitly declared proper material frame. The directions must share the original caloric identity, pure phase, temperature-axis convention and all fixed source context. The adapter uses the intersection of their original temperature and energy domains, preserves original source lines and receipts, and refuses extrapolation. Preparation compares fixed context once with bounded cancellation; evaluating a point or interval performs no material selection, table copy or source re-resolution.

Every crossed original conductivity knot contributes to each directional interval mean before rotation. The canonical `ConductivityCurve::mean_conductivity_w_m_k` accumulates positive segment means with convex width weighting, so representable means survive overflow or underflow of the dimensional Kirchhoff integral. Equal endpoints use the original point law; reversed intervals use identical arithmetic.

With `enthalpy-coupling`, the same adapter implements the original caloric storage trait. Its temperature, density and entropy queries retain the narrower source energy support. The binding original source owns each limiting energy endpoint, preserving its admitted one-sided temperature tangent at melting onset. No apparent heat capacity, new caloric identity, phase-mixture conductivity or inferred orientation is introduced. The adapter composes with the existing steady, sensible and caloric tensor drivers and original interface contacts.

Two actual production conductivity-mean tests passed under the pinned Rust compiler, including extreme representability and crossed-knot/support cases. Added adapter regressions cover source intersections, explicit rotation, original-context and chart admission, cancellation, and the pure-solid caloric endpoint; these integration definitions were parsed but not executed in the dependency-constrained workspace.

## Adaptive original-caloric tensor time integration

The original caloric tensor owner admits bounded adaptive SDIRK2(1) windows through `advance_caloric_tensor_conduction_adaptive`. Every completed attempt uses the same two native spatial stages; the embedded first-order energy discrepancy controls dt with the existing square-root thermal controller. Its scale uses changes in original h/u, including latent plateaus, never absolute energy datums. The embedded coordinate is never passed to the material inverse. Only the native second-order endpoint and native accepted stage heat are physical; the discrepancy is a local estimator, not a certified endpoint, global or channel-quadrature bound.

Forcing is sampled at actual native stage times. The callback is deterministic and read-only with respect to physical history and may borrow constant setup/prescription arrays. Strictly increasing declared transitions bound every attempted interval and mark endpoint samples as explicit left limits. Prescribed masks remain constant across all attempts in each smooth segment. Minimum/maximum dt, a maximum attempt count, and the native shared Newton/Krylov/backtracking caps bound work; only temporal rejection and native SolverNotConverged permit retries. Invalid source data, cancellation and other native refusals are fatal.

Rejected state and heat never enter retained physical records. Final stored energy is independently recomputed from initial/final original energies and fixed reference masses. Signed net accepted heat is compensated independently of the displayed source and outflow totals; closure subtracts the compensated storage and net-heat accounts before rounding away either correction. Both net whole-window closure and the sum of absolute accepted native closures must meet the single configured joule allowance. Any late error or cancellation discards the private accepted prefix and leaves caller inputs unchanged.

Five tests of the actual production embedded-error and heat-account arithmetic passed under the pinned Rust compiler. Five added native integration definitions cover diffusion refinement, latent energy and datum shifts, transition left limits and thermostat changes, bounded attempts and minimum dt, and fatal source/mask/cancellation outcomes. These integration definitions and public exports were parsed; full native integration execution was unavailable in the dependency-constrained workspace.

### Full-tensor diffusion on native open 3D flow

`OpenBoundaryGrid3::step_serial_with_thermal_tensor_diffusion` transports heat
with the original prepared `DiffusivityTensor<3>` and borrowed uniform or
per-cell `DiffusivityField`. It uses the existing D3Q7 flux block
`S = (I/2 + alpha_lattice/(1/4))^-1`, retaining off-diagonal transport in all
three Cartesian planes. Physical input scales as `K/(rho*Cp) * dt/dx^2`
with the carrier's fixed, common reference volumetric heat capacity. The
carrier's scalar tau continues to relax even kinetic moments. No tensor is
installed persistently or substituted into later scalar calls.

The original source half-kicks, advective-flux history, native force-corrected
velocities, actual wall/reservoir heat, finite-wall storage and final joint
publication remain shared with scalar transport. Solid grid cells carry no
thermal population inventory; explicitly configured immersed bodies retain
separate shared storage, while unassigned obstacle links retain the original
boundary law. Per-cell field dimensions are admitted before any candidate staging.
Refused or cancelled attempts publish no heat, momentum, wall state or time.
The original scalar collision path retains its arithmetic. G0 tests compare
the isotropic limit for D3Q19/D3Q27 flow with reversed reservoirs and sources,
and exercise late rejection/cancellation plus unchanged scalar retry. G1 tests
resolve signed oblique Fourier decay and spatial refinement in each Cartesian
plane against `exp(-t k^T alpha k)`.

This is anisotropic temperature diffusion, not variable-capacity conjugate
enthalpy or an inferred contact-resistance law. SPD alone establishes neither
population positivity nor full coupled-flow stability. Native wall and outlet
approximations retain their original scope; spatial/time refinement remains
the caller's responsibility.

### Tensor diffusion with native Boussinesq feedback

`BoussinesqOpenGrid3::step_serial_with_tensor_diffusion` carries the same
borrowed full thermal tensor into the original native temperature/force/flow
transaction. The transported temperature, including the final source half-kick,
drives the actual endpoint Boussinesq force. Off-diagonal conduction therefore
changes buoyant motion through the existing force owner and consistency audit.

Accepted and candidate force epochs, trapezoidal body impulse, reservoir heat,
finite-wall storage and wall conduction retain their existing accounting.
Temperature, momentum, walls, force history and time publish jointly; a rejected
or cancelled attempt retains the complete accepted state. The scalar entry
point retains its original arithmetic and no tensor persists into later calls.

Three regression definitions compare the isotropic limit on both native
momentum lattices, verify anisotropy changes endpoint force while closing the
fluid-plus-wall heat balance, and cancel at the final publication checkpoint
before an exact retry. This keeps the existing opening-velocity split and
small-density-departure Boussinesq assumptions. Frozen tensor coefficients and
a common reference heat capacity do not supply variable-storage enthalpy,
total-energy evolution, or a new temporal-order guarantee.

### Binary molecular diffusion on accepted carrier mass

`MaterialAdvectionStep::diffuse_binary_species` advances two original
constituent inventories on the SAME accepted carrier mass. It can consume
`advect_species::<2>`'s endpoint or a previous diffusion endpoint. Prescribed
per-cell molecular diffusivity is held fixed; it is never inferred from
viscosity or from a species index. The mass-frame law is
`j_A = -rho D grad(Y_A)`, `j_B = -j_A`, with density from accepted mass
and native voxel volume.

Every conducting fluid-fluid face uses the two half-cell resistances in
series, including graded Cartesian widths and spatial jumps in `rho D`.
Stagnant faces still diffuse. Exterior, solid and zero-diffusivity interfaces
are impermeable. Uniform mass fractions are stationary even when density
varies. The domain must be the original geometry: dimensions and fluid mask
are checked, while identical physical widths remain the caller's obligation.

Each pair's exponential relaxation approaches its mass-weighted equilibrium.
Forward/reverse half sweeps form a symmetric second-order diffusion splitting
for frozen coefficients. Both original constituent inventories receive
opposite transfers; neither clipping nor complementary-species reconstruction
is used. The less abundant constituent supplies the gradient to preserve
trace transport when majority fractions round to one. Prepared face factors
keep exponentials outside repeated sweeps. Work is
`O(cells + faces*substeps)`; storage is `O(cells + faces)`.

Independent constituent-mass and carrier-partition budgets admit the complete
result. Work/storage limits and bounded cancellation polling precede publication;
carrier, inputs and accepted state remain unchanged on refusal. Eight regression
definitions cover exact unequal-density pair relaxation, uniform fractions,
barriers/solids, tiny traces and species permutation, temporal/spatial convergence,
strong-step positivity/conservation, and invalid-data/budget/cancellation refusal.

This is an isothermal, fixed-carrier binary Fick component. It supplies no EOS,
partial-enthalpy transport, nonideal thermodynamic factor, Soret/barodiffusion or
arbitrary multicomponent closure. Large stable sweeps still require accuracy
refinement. The original donor advection followed by diffusion remains an
operator split with first-order overall accuracy; second order refers only
to the frozen diffusion sweep composition.

`MaterialAdvectionStep::diffuse_species_adaptive::<N>` adds automatic uniform
temporal refinement of that original operator. It compares complete n- and
2n-substep solutions from the same original inventories, accepts the finer
positive conservative state, and reuses a rejected fine state as the next
coarse comparison. Each constituent has its own absolute mass-fraction
tolerance; relative scaling uses the local initial/fine concentration, so
trace species do not inherit a bulk absolute allowance. The per-cell estimate
is `|Y_fine - Y_coarse| / 3`, from the second-order asymptotic regime. This is
not a certified global error bound and does not measure spatial, material or
advection-splitting errors. No Richardson extrapolation alters physical masses.
Aggregate work includes every discarded resolution; storage admission includes
the live coarse comparison alongside the original solver. Refusal or
cancellation returns no partial interval. Independent original constituent and
carrier-partition budgets remain mandatory at every resolution.

### Conservative multicomponent turbulent mixing

`MaterialAdvectionStep::diffuse_species::<N>` extends the same accepted-carrier
diffusion owner to 1..=128 constituents sharing a common diffusivity. Each
original constituent inventory evolves independently under
`j_k = -rho D grad(Y_k)`. One prepared face conductance and exponential factor
serve every constituent. The original binary entry and trace-safe complementary
arithmetic remain available. All components retain separate mass-defect budgets;
the sum is audited against the immutable accepted carrier, with no clipping,
normalization or reconstructed remainder species.

`FvFlow::diffuse_turbulent_species` consumes the actual momentum eddy viscosity
from the original LVEL or transient SGS flow, using an explicit positive
turbulent Schmidt number: `D_t = nu_t / Sc_t`. This is the gradient-diffusion
mass-flux closure in NASA/TM-2022-0008776, equation (46). The thermal Prandtl
number and LVEL heat-wall secant do not select species transport. Source-resolved
molecular viscosity remains with the momentum owner. The original accepted
mass divided by native cell volume supplies density, including graded grids and
spatial density jumps. Coefficients are frozen for the requested mixing interval.

`FvFlow::diffuse_turbulent_species_adaptive` connects the native momentum
eddy field to the original conservative temporal-refinement owner.
`FvFlow::diffuse_molecular_and_turbulent_species_adaptive` additionally consumes
one explicit nonnegative common molecular diffusivity per cell. The law is
`j_k = -rho (D_molecular + nu_t/Sc_t) grad(Y_k)`: molecular transport survives
zero SGS viscosity. The two cellwise contributions are added BEFORE the
original graded half-cell resistance is formed; adding independently
harmonized molecular and turbulent face coefficients is a different interface
model. Molecular D is never inferred from viscosity or the thermal Prandtl
number. Both coefficients must vanish in solids; a zero total coefficient
blocks the face. Negative/nonfinite inputs and unrepresentable positive
transport refuse without clipping or an invented diffusion floor.

Both adaptive entries retain per-species concentration tolerances, all
coarse/fine work in the aggregate budget, and the original mass/partition
accounts. Storage includes the combined D field and the simultaneously live
coarse comparison. The finer conservative state is authoritative; no
Richardson extrapolation alters its inventories. The temporal estimate is
asymptotic and does not certify spatial or model error. The molecular entry
admits binary Fick diffusion or an independently justified COMMON-D molecular
model; it does not stand in for unequal-species or Maxwell–Stefan transport.

The common face map preserves each species' global mass and the cell carrier
partition to the caller's declared budgets. Uniform mass fractions are
stationary at variable density. Zero eddy viscosity is an exact identity for
the turbulent-only entries; molecular diffusion can remain active.
The default entries keep exterior and solid faces impermeable; explicitly
authored surface reservoirs use the boundary entries below. Advective
inlet/backflow composition remains with `advect_species`. Work is
`O(N * (cells + faces * substeps))`, storage `O(N * cells + faces)`;
the checked component-work count and temporary diffusivity buffer are included
in the work/storage account. Cancellation or failed admission publishes no
partial inventory and does not modify the carrier or flow.

G0/G1/G3 regressions cover a four-species analytic pair on graded unequal-density
cells, trace transport and species permutation, the common-diffusivity Neumann
mode's second-order temporal refinement, null controls and resource refusals.
The native transient Smagorinsky-to-three-species test exercises the actual
momentum closure and accepted material flow. Separate coupling checks establish
Schmidt scaling, independence from the heat-wall policy and refusal atomicity.

This is an isothermal common-diffusivity mixing component. It does not supply
unequal molecular diffusion, partial-enthalpy transport, reaction-rate closure,
scalar variance or a species wall function. Matching original physical grid
widths, flow mask and the flow epoch remains the caller's responsibility. The
original turbulence model's applicability and resolution limits remain in force.


### Reservoir-driven species diffusion

`MaterialAdvectionStep::diffuse_species_with_boundaries` adds explicitly
authored fixed-composition reservoirs at exterior or fluid-solid faces.
Each face supplies unit-sum mass fractions and a common finite nonnegative
film resistance in m² s/kg. With the original accepted density and native
half-cell width, `G = A / (w/(2 rho D) + R_film)` in kg/s. Zero film resistance
prescribes composition at the actual face; zero cell D still blocks exchange.
Boundary identities must be unique and ordered by `(cell, face)`; interior
fluid-fluid faces cannot acquire a second transport owner. Unlisted faces
remain impermeable. A solid-facing reservoir does not evolve solid storage.

Exact positive cell/reservoir relaxation maps surround the original interior
half sweeps in canonical forward/reverse order. Signed, gross inward and gross
outward constituent transfers are retained at each face. Original constituent
mass admission subtracts the actual signed external transfers; every boundary
map's nonzero summed carrier transfer is charged to the partition budget.
Both boundary and interior updates count toward work, and prepared coefficients,
transfer accumulators and face readbacks count toward storage.

`diffuse_species_with_boundaries_adaptive` reuses the original complete-interval
coarse/fine refinement. Besides per-cell concentration tolerances, each species
has an explicit integrated-transfer allowance in kg. The estimate sums, over
faces, the largest signed/inward/outward transfer difference divided by three;
opposing reservoirs cannot cancel unresolved transport. The retained fine
inventories and transfer accounts are authoritative, without extrapolation.
These are asymptotic temporal estimates, not certified spatial/model bounds.
This is an isothermal common-D reservoir law, with no inferred selective
membrane, sorption equilibrium, finite reservoir/solid storage or enthalpy flux.

### Adaptive accuracy for caloric reservoir diffusion

`diffuse_species_with_enthalpy_and_boundaries_adaptive` refines the original
fixed-coefficient reservoir caloric transaction from the same initial stores.
Every resolution retains the same source-owned face laws and constituent order.
The original concentration criterion is joined by independent per-cell specific
enthalpy and per-face species/energy integral criteria. Cell-energy differences
use retained increments divided by the original carrier mass; face differences
use retained signed and both gross-sign accounts. Their largest absolute
difference is summed across faces and divided by three. Large formation data,
opposing reservoirs and rounded scalar readbacks cannot erase those estimates.

The accepted fine inventories and original receipts are unchanged. All trials
share the cumulative pair-work budget; retained coarse fields and original
per-cell/per-face accounts enter the storage allowance. These are asymptotic
temporal estimates of the frozen model, not spatial, source-refresh or phase
transition guarantees. The exact compensated-comparison helper passed its
focused native test. Six analytic/refinement/refusal regression groups are
authored but unexecuted because the dependency build is blocked before fs-lbm.

### Native open 3D local SGS momentum and heat

`lattice3::open_boundaries::sgs` prepares Smagorinsky, WALE or Vreman
on the original D3Q19/D3Q27 reservoir geometry. The accepted, force-centered
velocity field uses the native full gradient: quadratic one-sided velocity
reservoir slopes (linear when the second donor is solid), the pressure
reservoir's zero-normal-gradient donor convention, and unchanged interior,
periodic and stationary-wall stencils. Open faces are never interpreted as
impermeable halfway walls. Shared local SGS admission and collision-family,
bulk/higher-mode preservation are retained. The original molecular collision
or a fresh explicit per-cell molecular field is the baseline on each call;
previously augmented collisions must not be reused as molecular transport.

`OpenSgsWorkspace3` retains the actual prepared eddy viscosity rather than
recovering it by subtracting rounded collision viscosities. Its explicit
numerical payload budget covers five f64 values and one native collision per
cell, with actual Vec capacities checked. Serial and pooled native momentum
steps consume these prepared collisions without adding persistent closure
state. Preparation performs no adapter-owned allocation.

`step_serial_with_thermal_sgs` couples that same local eddy viscosity to the
original native open-flow/D3Q7 heat/finite-wall transaction. Its explicit
positive turbulent Prandtl number gives alpha_t = nu_t/Pr_t. Since the actual
D3Q7 scalar sound-speed squared is 1/4, it adds 4*alpha_t to each supplied
molecular thermal tau. The caller supplies a fresh complete molecular thermal
tau field. A separate nonnegative finite alpha_t/alpha_molecular maximum is
an admission limit, never clipping. Nonfinite or underflowed positive modeled
transport refuses; ordinary finite floating-point addition may round a tiny
positive contribution back to the molecular tau. Exact zero SGS preserves
both molecular collision and thermal-tau bits. The coupled workspace adds
one f64 per cell to the momentum workspace.

The existing native owner retains reservoir reversal, scalar-population heat
exchange, symmetric source kicks, advective history, finite-wall storage and
wall conduction. All coefficient admissions precede staging, and its final
heat balance and cancellation gates precede joint publication. A refusal may
change scratch only: accepted momentum/thermal populations, finite walls,
receipts and step counters remain unchanged. Modeled mechanical transfer
remains diagnostic, with no invented viscous-heating source. This is passive
fixed-reference-capacity low-Mach temperature transport with frozen incoming
coefficients, not compressible enthalpy, a dynamic exterior test filter,
nonreflecting-outlet physics or a qualified turbulent-jet result.

Four native regression groups cover D3Q19 TRT and D3Q27 cumulant/full-gradient
reservoir preparation for all three models; exact zero-model joint trajectories
with sources and conducting finite walls; active native composition and an
independently measured fluid/solid heat inventory; and warmed-state budget,
ratio, late thermal, cancellation and retry behavior including positive
sub-ulp modeled transport. These regressions are authored and parsed/formatted;
the native integration suite was not compiled or executed in the constrained
implementation workspace.

## Reciprocal momentum on evolving material liquids

The sealed `MaterialFreeSurfaceMomentum3` interface binds the existing
`MaterialFreeSurface3`, `MaterialFreeSurfaceRheology3` and
`MaterialFreeSurfaceEnthalpy3` owners to the same compact native impulse
transaction. It exposes the original reference-density, length, velocity and
scenario scaling together with a borrowed native wet-cell stage. Implementations
delegate to the original physical speed/unit admission and retained liquid
inventory; composing solvers cannot substitute another material or flow owner.

This instantaneous mechanical source preserves canonical density, actual partial
mass, held force frame and physical nonequilibrium stress. The thermal and
rheological owners retain their accepted source samples, caloric/phase/rheology
history and chronology. No collision, constitutive query, heat deposition,
geometry update or clock advance is inferred. Dropping a refused stage publishes
nothing, and no complete material owner or population field is cloned.

The interface enables reciprocal finite-body force/couple and dissipation
checks in `fs-couple` before native publication, with the same source binding
for frozen, non-Newtonian and thermally evolving liquids.

## Tangential moving walls in native open flow

`OpenBoundaryGrid3` now admits prescribed tangential exterior-wall velocities
for both D3Q19 and D3Q27. The existing native moving reflection and compensated
relative-velocity impulse/work reduction serve closed and open grids. Reservoir
donors and reservoir exchange use the same moving-wall candidate. At open/wall
seams, reflected links stay wall owned; the reservoir cannot conceal wall mass
injection. Every reflected set passes a local zero-mass-source check before
admission. Streamwise moving-wall through-flow is supported, including its seams.
Normal motion and incompatible velocities at intersecting wall faces refuse;
the fixed geometry does not represent swept volume or fresh cells.

Boundary rebinding and serial/pooled steps accept wall-speed updates without
resetting populations. `last_wall_exchange` measures all exterior wall links
after an accepted moving interval, separately from reservoir transport. It
retains wall, population and moving-mass impulses and their angular counterparts,
plus prescribed wall work. Refusal preserves the accepted populations, targets,
time and last receipt. Spatial-force and thermal publication retain this same
transaction. Rheology and SGS gradients sample the actual wall velocities.
The thermal carrier supports stateless normal wall laws on tangentially moving
faces; finite moving-wall storage refuses because solid heat advection is not
represented. Mechanical wall work is not silently converted to heat.

Focused tests cover uniform co-motion with through-flow in all orientations,
separate wall/reservoir momentum closure and seam populations, physical wall
shear, force-epoch publication, and native-kernel parity through analytic
Couette convergence. The linear moving reflection retains the native low-Mach
model; no new high-Mach or moving-geometry accuracy claim is added.

## Spatial velocity traces on native open flow

`OpenBoundaryGrid3` and `ForcedOpenBoundaryGrid3` accept spatial velocity
profiles on all six native face orientations for D3Q19 and D3Q27. They reuse
the specialized owner's `InletFlux3`, profile report and refusal types. Callers
may resample at successive physical times without resetting populations or
force history; two face-sized vector buffers are reused, with no profile
callback or allocation during stepping. Sampling/refusal leaves accepted
populations and time unchanged. A successfully installed control remains held
if a later fluid step refuses.

Optional outward-normal mean control uses compensated deterministic reduction.
Both sampled and corrected vectors respect the native grid/collision Mach
envelope; an optional inclusive speed limit applies to the final corrected
trace. Mixed wall/open seams retain their original per-link ownership and
refuse mean-normal control because their complete velocity is not imposed.
Unrelated wall/reservoir updates preserve profiles. An attached profile must
be cleared before changing its uniform fallback; clearing restores that held
fallback without changing physical state.

The original force-corrected donor stress and regularized reconstruction consume
the held vectors directly. Momentum receipts measure the actual operator, while
thermal advection, rheology and SGS consume the accepted physical velocities.
No separate flow or turbulence source model is introduced. Focused cases cover
time updates across orientations/collision families, momentum and mass closure,
uniform-profile bit parity, native Mach/final-speed refusal, late cancellation,
moving-wall seam ownership, and spatial force epochs. The requested mean is a
volume-velocity control, not an assertion of exact mass flux.

## Native open-grid spatial force epochs

`lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3` binds an explicit per-cell force density to the original native D3Q19/D3Q27 open-grid populations. Initial attachment changes only the native half-source coordinates: density, physical velocity, nonequilibrium stress, accepted time and reservoir/pressure controls survive. There is no public inner-grid escape with stale uniform-force macroscopic readback.

A step collides with the accepted opening force and interprets streamed donors and targets with their own endpoint force. Donor physical stress uses the donor source; reconstruction attaches the target source. The fluid body impulse is the trapezoidal sum of the two force epochs, separate from measured reservoir and wall transport. Ordinary steps do not perform an extra source-coordinate reframe. Both source fields must satisfy original collision capability, finite-value, shape and solid-zero admission. Cancellation, Mach or reconstruction failure leaves accepted populations, force history and time unchanged. The wrapper retains original native collision kernels and preallocated grid scratch; holding a spatial force allocates no per-step adapter storage.

Focused checks execute the actual generic-poll native owner cores across five D3Q19 and six D3Q27 collision configurations: physical velocity/stress attachment, periodic trapezoidal impulse, bitwise constant-field equivalence, all three open axes with distinct donor/target sources, measured reservoir closure, and late refusal/exact retry. These extracted-production checks exclude the thin `Cx`/pool wrappers and do not substitute numerical kernels. They do not constitute a full Cargo integration run or a new finite-force accuracy claim for the native collision families.

## Native open-flow Boussinesq temperature feedback

`lattice3::open_boundaries::thermal::buoyancy::BoussinesqOpenGrid3` owns native D3Q19/D3Q27 momentum, D3Q7 temperature, finite thermal walls/conduction and both force epochs as one accepted state. Its explicit affine law is `F_b = -rho_ref * beta_theta * (theta - theta_ref) * gravity_lu`, plus an authored per-cell mechanical background that is held or supplied at the next accepted time. Here `beta_theta = beta_SI * delta_T` and `gravity_lu = gravity_SI * dt^2 / dx`; the law uses fixed reference density while native macroscopic velocity retains the actual population density. Attachment preserves physical momentum/stress and existing thermal history. This is fixed-reference-capacity heat transport with Boussinesq momentum feedback, not variable-density mass transport or a total-energy equation.

Each step transports heat at the accepted force-corrected velocity, previews the endpoint force from transported temperature plus the remaining source half-kick, and performs one original native momentum step. Collision uses the retained opening force; donor/target reconstruction and endpoint readback use their local candidate force. The final scalar source kick uses the admitted native endpoint velocity. Its actual temperature must reproduce the applied force within the caller's absolute component tolerance; the accepted representation retains that applied force. Scalar splitting retains its original temporal-order limits. The separately reported trapezoidal body/buoyancy impulse is not deposited as sensible heat.

The declared relative-density-departure envelope admits cell endpoints/source previews, authored reservoir/fixed-temperature wall data and actual candidate fluid traces at reservoirs and external/internal walls. Reservoir receipts retain the extrema of their reconstructed incoming/outgoing population pairs, including outflow/backflow. Finite-solid bulk temperatures remain governed by the thermal model, not the fluid density law. A refusal or cancellation publishes none of the flow, force, thermal history, wall state, time or receipts. Three explicitly budgeted force fields and the existing native/scalar/solid scratch are reused without per-step adapter allocation; no iterative heat re-solve is needed.

Focused checks execute five actual coupled-owner cases: analytic uniformly heated boxes on both native lattices with reference density distinct from population density; thermal reversal with zero net buoyancy impulse; a conducting-wall open channel with measured two-way feedback and independent fluid-plus-solid thermal closure; final-poll/domain refusal with bitwise exact retry; and fluid-boundary domain enforcement despite admissible mixed cells. Two direct law tests cover signs, temperature-datum shifts, anomalous expansion, zero-expansion background preservation and explicit domain refusal. The numerical-owner tests use extracted production cores and real native/math/sparse kernels, excluding only execution-context wrappers. Full Cargo/public-context integration and the planned developing-channel Nusselt/refinement reference remain unverified by these tests.

### Mixed buoyant and prescribed-force convection

`BoussinesqOpenGrid3::advance_serial` and `advance_serial_with_tensor_diffusion`
accept a caller-owned `background_at_end: &mut Vec<[f64; 3]>`. The field is the
mechanical force density at the next accepted lattice time, before adding
thermal buoyancy. It covers every original cell, is finite and vanishes on
solid cells. Spatially nonuniform and time-varying backgrounds use the same
native D3Q19/D3Q27 momentum and D3Q7 heat owners as held-background steps.

Collision uses the retained total force `F_n = B_n + F_b(theta_n)`.
Reservoir reconstruction, endpoint physical velocity and the final scalar
source half-kick use `F_(n+1) = B_(n+1) + F_b(theta_transport + Q/2)`.
The original final force-consistency gate checks the actual accepted
temperature against this endpoint background. Thermal advection retains its
opening physical velocity and existing source split; no additional temporal
order is claimed.

`BoussinesqOpenStep3::background_impulse` reports the fluid-cell sum of
`(B_n + B_(n+1))/2`. `body_impulse` uses the corresponding total-force
endpoints, and `buoyancy_impulse` subtracts the mechanical background using
compensated accumulation. Wall and reservoir momentum remain separate.
Mechanical and buoyant work are not deposited as heat; the original reservoir,
finite-wall, immersed-solid and thermal-source accounts remain authoritative.

After the final cancellation gate, accepted flow, total force, thermal state,
wall histories and mechanical background publish together. Success exchanges
the caller's endpoint buffer with the old accepted background for reuse;
refusal or cancellation preserves both the input buffer and accepted epoch.
The owner retains its existing three force fields without a per-step force
allocation. Caller endpoint storage is outside `maximum_force_storage_bytes`,
which bounds retained field payload. An unchanged endpoint background retains
held-background arithmetic. Focused regression definitions cover affine
forcing/heating on both lattices, identical-endpoint equivalence, spatial
forcing with tensor/finite-wall heat closure, and refusal with exact retry.

## Constituent mass on accepted material transfers

`MaterialAdvectionStep::advect_species` transports 1..=128 fixed-order
constituent inventories on the retained native three-dimensional material
advection routes. Per-cell constituent **mass** is authoritative; its initial
partition is checked against the accepted step's original mass and its endpoint
partition against that same step's represented endpoint mass. Incoming boundary
planes and local sources have separately authored mass fractions selected by
each retained route's actual origin. Extraction and backflow use their original
donor/origin rules. No second density solve, inferred partial volume, or
composition-dependent EOS is introduced.

Each species reuses the original row-scaled mass-times-specific-quantity donor
matrix and the existing BiCGStab recurrence. One native ILU(0) factor is reused
across species, with only one scalar solve workspace live at a time. The caller
bounds component work and payload storage. Per connected component, fraction
amplitude scaling keeps an entirely trace species representable in the native
Krylov norm; absent species retain the original exact zero-RHS path. This does
not promise arbitrary relative accuracy for a tiny local tail of a bulk species.

Paired accepted transfers determine endpoint constituent masses. Independent
readback substitutes those actual masses into the unscaled extensive donor
equations, checks each species' global boundary/source balance, and audits the
initial, endpoint and route partitions against the accepted carrier. Each
species has its own absolute mass-defect budget, including trace species;
partition errors have a separate declared kilogram allowance. Pure-material
derived fractions can differ from one by admitted carrier/constituent roundoff.
There is no clipping, seeded trace or composition renormalization. Negative or
nonrepresentable inventories/transfers and exceeded budgets return no candidate.
Cancellation, including after completed components, leaves all caller data and
the accepted carrier step unchanged.

This is the native first-order implicit donor advection scope. Species diffusion,
reactions, mixing heat, thermodynamic composition feedback and a higher-order
advection scheme remain with their respective future/material owners.

## Thermal evolution coupled to resolved Carreau film momentum

With `implicit-enthalpy-coupling`,
`film::transient::thermal::ThermalCarreauFilm` binds the retained original
wall-normal velocity to the original identity-bearing caloric timeline.
Ordered fluid cells must form the same straight, uniform wall-normal mesh:
each cell retains volume `A*dy`, mass `rho*A*dy`, the same original caloric
material and the same accepted time. Existing solid or wall cells may share
the native thermal timeline. No capacity, conductivity or liquid phase
occupancy is inferred from a viscosity law.

Each native backward-Euler thermal residual, derivative probe and independent
final evaluation restarts the SAME accepted momentum profile. Interior shear
face temperatures interpolate the original native caloric cell temperatures;
the callback supplies the actual wall and free-surface fluid-side traces
from its original interface/conduction model. A solid bulk temperature is not
automatically an interface temperature. Those endpoint source temperatures
are held during the original midpoint momentum solve. This composition is
first order in coupled time, even though momentum at a prescribed fixed
temperature retains its midpoint scheme.

The original local viscous-loss partition supplies each fluid cell's heat
exactly once. Signed body, moving-wall and free-surface work remain mechanical
exchange; neither kinetic-energy change nor algebraic defects are deposited
as heat. Original native thermal acceptance precedes a separate combined
`delta(thermal + kinetic) - other heat - mechanical work` gate. Other thermal
power is measured directly by the native owner before viscous heat is added;
the combined account uses the independently measured stored-energy change.
Reports retain actual deposited heat, its partition residual, kinetic change,
mechanical work and thermal/nonlinear residuals separately.

The owner reuses its momentum trial arrays, temperature buffer and local heat
buffer. Every thermal trial resets only the original accepted velocity/time;
rejected Newton candidates do not become physical histories. Endpoint and
aggregate momentum constitutive budgets bound nonlinear work. Native thermal
and Krylov allocations remain. Cancellation is latched for the whole advance.
No accepted velocity, thermal state or time changes unless all final gates
pass; an explicit `into_parts` handoff permits native events and later rebinding.

This is a locally parallel, fixed-thickness, fixed-reference-density liquid
film. The thermal source must remain in its admitted liquid domain.
Streamwise transport, evolving interfaces, thermal expansion and phase
occupancy are not introduced by this coupling. Nonlinear convergence alone
does not establish temporal or spatial accuracy.

The coupling's default Newton controls use an inner forcing floor of `1e-4`,
appropriate for its finite-difference Jacobian and caloric/constitutive roundoff.
The original outer residual and physical energy gates remain independent;
explicit caller-supplied solver controls are unchanged.

All 37 focused source-linked film tests passed, including four coupled
regressions for viscous self-heating and original-source feedback, an independent
two-cell conduction/bath balance, first-order refinement against independent
coupled momentum/thermal integration, and joint refusal/cancellation/retry.
Executed code includes the original material laws, math, tridiagonal, thermal
backward-Euler, finite-difference and Newton/Krylov kernels. The harness omits
unrelated crate modules and uses an uninhabitable context type for unexercised
runtime entrypoints. These results do not establish full-crate Cargo or
asupersync runtime integration.

## Resolved transient Carreau films

`set_face_temperatures` admits an explicitly supplied wall-normal temperature
field at all N+1 original shear faces. Every face retains the SAME thermal
Carreau source and physical rate ceiling, with its own prepared temperature
shift. Interior stresses/tangents and the free-surface stress inverse use those
actual local preparations; no mean temperature or fitted apparent viscosity is
substituted. Uniform fields retain the original compact path. Spatial updates
reuse two prepared buffers after first use and publish only after complete
source/scaling admission and the final cancellation check. Updating temperature
never replaces velocity history, changes density/thickness, or advances time.

The Newtonian, original-source-bound Newtonian and Carreau `step_with_heat`
entries share one shear-segment partition. Each exports a nonnegative
viscous-energy increment per original
wall-normal cell, in J/m² of wetted wall area. Interior shear-segment work is
split equally between its two adjacent cells; the wall and surface half-segment
losses enter only their adjacent cells. The sum is the original momentum work
account's viscous heat, within summation roundoff. The surface contribution uses
the imposed traction and retains the separate constitutive-inverse defect.
The caller supplies its actual represented wall area and original caloric
owner; neither a heat capacity nor a temperature rise is inferred. These
increments already exclude kinetic storage, so body, wall and surface work
must not be deposited again. Heat output, velocity and time publish together;
shape, source, solve, speed and cancellation refusals retain every accepted
value. Reusable O(N) scratch avoids momentum-step workspace allocation.

G1/G3/G4 spatial-field tests compare against independent original-law stress
bisection and continuum integration under a linear temperature gradient and
opposing surface traction; they check spatial refinement, uniform-path parity,
local heat partition and global work, original face/source retention, thermal
rebind refusal, and cancellation/retry including the caller's heat buffer.
Native Newtonian checks also recover exact signed Couette cell heating and
bit-identical momentum/report behavior with and without heat export.
The temperature field is prescribed over each momentum interval. These APIs
provide spatial source feedback and conservative heat transfers, not a second
thermal inventory, conductivity closure or automatic thermal time integrator.

`film::transient::TransientCarreauFilm` retains the full wall-normal velocity
history under signed body acceleration, translating-wall velocity and imposed
surface traction. Each implicit-midpoint face evaluates the original thermal
Carreau stress and positive stress tangent. The shared `fs-material`
`PreparedCarreau` owner retains the original temperature shift and constitutive
inverse; no fitted apparent viscosity or summed steady profiles are substituted.
Bounded Newton/backtracking uses the native tridiagonal solve and reusable O(n)
arrays. Accepted trials retain their constitutive/Jacobian evaluation.

Newtonian and nonlinear films share one finite-volume momentum and work account,
including both boundary half-cells. The nonlinear surface inverse exposes its
stress defect separately. Initial, midpoint and final rates remain within the
authored source-supported ceiling; velocity, source, iteration, work and final
cancellation refusals leave accepted history and time unchanged. A prescribed
temperature update retains momentum, density, geometry and the physical rate
ceiling while re-preparing the original source and diffusion scale.

G1/G3/G4 tests cover the original Newtonian limit, independent constitutive
integration and spatial refinement, nonlinear time refinement, reversing and
moving-wall work, exact rest/Couette limits, temperature continuation and atomic
retry. This is fixed-thickness, locally parallel generalized-Newtonian flow.
It adds no elastic memory, thermal evolution, free-surface geometry, entrance
flow, capillary waves or turbulent-transition prediction. Midpoint time and
wall-normal resolution remain accuracy requirements.

## Resolved transient laminar films

`film::transient::TransientFilm` evolves the wall-normal velocity of the same
constant-property Newtonian film used by `NusseltFilm` and `ShearedFilm`.
The local equation is rho u_t = rho g + mu u_yy at fixed thickness, with
prescribed tangential no-slip wall speed and surface traction. Uniform finite
volumes and implicit midpoint retain startup, decay, reversal and counterflow;
drive changes do not replace the accepted profile by a steady parabola.

The original native tridiagonal solver uses O(n) work/storage, reusable scratch
and a cached factorization for unchanged durations. Every numerical traversal
polls cancellation at least every 128 cells. The caller supplies physical-speed
and original-equation residual limits. Invalid source/geometry/drive, failed
solve, exceeded limits or cancellation publishes neither velocity nor time.

Receipts separate body, wall and surface work, represented kinetic change,
nonnegative viscous heat, and momentum/work defects. The prescribed surface
stress reconstructs the actual boundary trace across its half-cell; that
half-cell's traction work and viscous loss occur exactly once. Heat remains a
reported transfer, not a second thermal inventory or inferred temperature rise.
`MaterialFilmLiquid::transient_film` retains the original liquid source,
viscosity basis, query and receipts beside this native history. No additional
material selection or apparent-viscosity fit is introduced.

G1 tests compare startup with a continuous mixed-boundary sine expansion,
measure second-order spatial and temporal refinement, and recover the driven
Nusselt/sheared steady limit. G0/G3/G4 checks cover signed work and momentum,
common-frame boosts, reversal, rest, domain refusal and exact cancellation
retry. These tests exercise the numerical owner; the material adapter also
has source-to-native parity regressions for both viscosity bases. Verification
of the numerical owner uses exact production film/math/tridiagonal sources;
full Cargo and source-adapter integration remain separate verification gates.

This is fixed-thickness, locally parallel, laterally wide Newtonian laminar
flow. It does not evolve a free surface, infer turbulent transition, resolve
entrance effects or change material temperature/composition. Midpoint may
oscillate on unresolved stiff modes; conservation and linear stability do not
replace temporal or spatial resolution checks for the forcing bandwidth.

## Source-resolved Carreau viscosity in nonlinear slender jets

`PeriodicSlenderJet::new_carreau` retains an original
`ThermalCarreauYasuda` law and fixed physical temperature, density and surface
tension. Each actual explicit or implicit RK stage queries the material at
`gamma=sqrt(3)*|du/dz|`, the invariant of incompressible uniaxial extension.
The conservative axial stress is `3*mu*du/dz`; the same sample supplies
`3*mu*(du/dz)^2` to the existing integrated viscous-loss account. No lattice
tau cap, fitted Newtonian substitute or second constitutive implementation enters.

The explicit viscous step bound retains both apparent and differential
viscosity. The apparent conductance is necessary near strongly thinning states:
a small tangent alone can allow an Euler stage to overshoot through zero strain.
Every explicit trial stage rechecks the bound and the original source support.
Source refusal or cancellation preserves accepted area, momentum, geometry and
time. The explicit path retains its original numerical workspace.

Regression definitions cover exact Newtonian-degeneracy trajectories, shifted
temperature, independent tensor stress and heat, conservative momentum/volume,
strong-thinning dissipation, the viscous timestep bound, and refusal/rollback.
These definitions require execution in the repository's Rust validation lane;
they are not a claim that a new convergence or performance campaign has run.

`ImplicitViscosity` retains the original IMEX-SSP3(4,3,3) tableau for Carreau
jets. Each fixed-area nonlinear stage solves the original conservative stress
with its exact positive differential tangent. The native cyclic LDL/rank-one
kernel resolves the symmetric positive definite Newton Jacobian in O(N) work.
Five additional arrays are prepared once for a non-Newtonian selection; successful
stepping allocates nothing. Exact Newtonian degeneracies retain the native direct
solve and arithmetic. `SlenderJetImplicitControls` bounds corrections, backtracks
and relative stage residual; the outer step budget bounds the number of stages.

The source-weighted mean is separated before nonlinear residual scaling, so a
uniform translation does not become the convergence tolerance for weak shear.
Backtracking re-queries the same original material; unsupported algebraic trial
rates may be rejected without clipping or extending source support. A search
with no supported candidate preserves its original source diagnostic. Actual RK
stage refusal, exhausted nonlinear work and cancellation roll back the full
advance. Returned momentum and source rates are rechecked after mean restoration.

Positive RK weights integrate physical viscous loss at the solved stage rates.
The fixed-stage kinetic balance distinguishes molecular loss from implicit
numerical damping; full IMEX endpoints retain a signed truncation remainder.
Third-order temporal behavior requires resolved dynamics and nonlinear residual
error below time-discretization error. The embedded second-order estimate remains
a local indicator, not an unconditional accuracy or monotonic-energy guarantee.
Focused regressions use an independent nonlinear alternating-face decay law,
original-stress/tangent differences, stage energy balance, large translations,
stiff shear-thinning evolution, source-domain backtracking and transactional retry.

Temperature does not evolve, loss is not deposited as heat, and there is no
elastic memory or independently calibrated extensional response. Original
slenderness/resolution limits and pre-pinch stopping remain.


## Prescribed axial temperature and local viscous heat in nonlinear jets

`set_carreau_cell_temperatures` binds one prescribed fluid temperature to each
original Eulerian axial cell. Cell i owns the extensional stress between its
two neighboring velocity faces. Every preparation retains the SAME original
thermal Carreau source, with a caller-authored positive physical invariant-rate
ceiling. The complete interval from rest to that ceiling must be supported at
each temperature. No inferred floor, shifted reference-support clipping or
alternative viscosity fit is introduced. The original uniform-temperature
path remains available, including a source supported only at rest.

The original `PreparedCarreau` owner evaluates each temperature shift once.
Every explicit stress, implicit Newton tangent, residual recheck and physical
loss then queries the original prepared curve at its actual local invariant.
A temperature prescription also checks the accepted strain field before
publication; failed source/allocation/cancellation leaves the preceding source
field, trajectory, policy, diagnostics and time intact. Only reusable solver
scratch may change. `prepared_carreau_cells` exposes the actual temperatures,
source law and physical rate ceiling; the original parameter getter retains
the original calibration and initial uniform temperature.

An exact Newtonian degeneracy uses the old direct constant-viscosity solve only
when all prepared zero-rate viscosities are equal. Nonuniform thermal viscosity
uses the same cyclic weighted Jacobian; a rate-linear law needs one correction.
Preparing a new nonuniform field after selecting implicit integration prepares
the missing nonlinear workspace once. No scalar mean temperature or mean
viscosity replaces a spatially varying original source.

`enable_cell_viscous_heat` prepares four reusable arrays for local stage power,
trial-step quadrature, whole-call candidate and accepted heat. Explicit SSPRK3
and IMEX use the SAME positive weights as their existing aggregate loss report.
Each cell accumulates `3*pi*dz*mu*A*(du/dz)^2` at actual solved stages.
Rejected attempts contribute no heat; a domain stop retains only completed
substeps. Error/cancellation preserves the preceding whole-call heat account;
a successful zero-duration call publishes zero. Successful steps allocate
nothing in this account. `cell_viscous_heat` reports the local joules, the
identical aggregate loss and their unmodified compensated partition residual.

Only physical viscous loss enters this account. Upwind dissipation, implicit
damping and algebraic defects are not deposited. Temperature is held until an
explicit prescription changes it; density and surface tension remain fixed.
There is no caloric evolution, thermal expansion or Marangoni stress. A future
thermal consumer must transport the changing cell mass and its advected energy
as well as this heat; depositing into a frozen-mass thermal mesh is not equivalent.

Focused regression definitions cover independent local stress/tangent and
power, uniform-field trajectory parity, variable-temperature Newtonian limits,
third-order nonlinear refinement of each local heat family, and transactional
source/work/cancellation refusal and retry. The independent alternating-flow
reference includes evolving areas and two original temperature-shifted laws.
These definitions are not a claim that native Cargo tests or a performance
benchmark have run.

## Conservative thermal evolution and viscous self-heating in slender jets

Under `enthalpy-coupling`, `jet::nonlinear::thermal::ThermalCarreauJet`
binds the native explicit Carreau jet to an original
`IsochoricHeatCapacityLaw`. The source's fixed density must match the flow;
its original nonlinear Cv(T), caloric support and exact segment integration
and inversion remain authoritative. An optional `IsochoricConductionLaw`
must retain that exact caloric identity.

One coupled SSPRK3 evolves native area, staggered momentum and extensive cell
internal energy. Eulerian mass is `rho*pi*A*dz`. Energy advection reconstructs
minmod-limited linear traces of the anchored specific energy and shares the SAME
signed native area-cell face flux q between its two cells. This includes
transport when the dual momentum flux happens to vanish. Constant specific
energy follows the identical mass update, including signed relative energies. Each actual stage
inverts its evolving caloric coordinate and prepares the original thermal
Carreau source at that temperature. The original invariant, extensional stress,
analytic tangent and physical viscous power remain the existing owners.

The transported coordinate is `m*integral(Cv dT, T_reference, T)`, with a
single fixed numerical reference temperature. No subtraction of large source
energy datums enters transport; changing the source's declared zero leaves
this trajectory unchanged. The original preserving-energy inverse retains
nonzero stored heat smaller than one temperature ULP without inventing a
temperature increment. Conductivity is queried directly at the already
resolved temperature through its original prepared source, avoiding a
redundant absolute-energy round trip.

Only local physical viscous power enters calorics, with the same positive RK
quadrature as the native aggregate and cell heat accounts. Optional axial
Fourier exchange uses the series resistances `dz/(2*pi*k*A)` of the actual
two half-cells. Shared face powers enter neighboring inventories with opposite
signs. The explicit thermal bound combines twice the gross outgoing mass rate with
the conductive row sum times the original global maximum dT/du. This sufficient
bound preserves mass-weighted convexity for transport with changing cell mass;
physical heating still requires admission by the original caloric source. Native advection,
capillary and apparent/differential-viscosity restrictions remain active.

Flow, energy inventory, prepared source field, time and cell heat publish
together. Source/rate, work, arithmetic, conservation or cancellation refusal
rolls back the entire call. A native geometry stop retains only completed
coupled substeps and restores the final admitted source profile. An explicit
absolute-joule allowance checks stored-energy change against physical heat.
Mass roundoff, the thermal residual and the combined mechanical/thermal
truncation residual are exposed separately; upwind loss is never deposited.
The explicit numerical arrays are prepared once, with O(N) reusable storage; the limited
traces require no additional array. Identical rounded stage temperatures reuse
their immutable prepared rheology. Caloric inversion traverses crossed original
source segments and is not claimed to have logarithmic cost for every
reference-relative query.

The explicit `advance` entry has third-order temporal accuracy for smooth admitted
semidiscrete solutions. Energy reconstruction is second order on smooth
monotone patches; limiting reduces spatial order near extrema. Density,
composition and surface tension stay fixed; the outer surface is insulated.
There is no radial temperature profile, expansion, phase change or Marangoni stress. Material incompressibility includes radial
strain, so pressure dilatation vanishes; native capillary work remains mechanical.
The original slenderness/resolution and pre-pinch limits still apply.

Eight jet regressions and two original-caloric/conductivity regressions define
coupled nonlinear-Cv evolution, temporal refinement, actual mass-flux energy
transport, uniform-temperature moving mass, source-datum invariance including
variable conductivity, source identity/density, geometry stopping, late
source/work/cancellation rollback, sub-ULP stored heat and original support.
An independently executed conservative alternating-mode reference gives
SSPRK3 error ratios about 0.1245 on timestep halving and thermal closure at
roundoff. A separate exact-cell-average periodic advection fixture covers both
flow directions and signed specific energies at 32, 64 and 128 cells. Its
independently computed L1 errors are 7.2683035e-4, 2.2051748e-4 and 6.0154975e-5
J/kg, with refinement ratios 0.3034 and 0.2728 approaching second order despite
limiter reduction at extrema. The corresponding donor error at 128 cells is
9.6669464e-4 J/kg. No new extrema or unaccounted heat appear in that independent
transport calculation. All changed Rust files parse with pinned rustfmt. Native
compilation and these Rust regressions remain unexecuted in the current
dependency/storage-constrained workspace; the independent oracle is not
native integration validation.

### Coupled implicit viscosity, conduction and self-heating

Under `implicit-enthalpy-coupling`, `ThermalCarreauJet::advance_implicit`
provides first-order IMEX Euler. Original limited conservative transport and
capillarity advance first; a jointly converged backward-Euler momentum and
internal-energy stage follows at the transported areas. Every residual resolves
the original relative caloric inverse, Carreau rate/temperature response and
direct-temperature conductivity. It does not freeze viscosity in the accepted
physical stage. An original held-temperature momentum solve supplies only an
algebraic stiff initial guess.

Analytic nonsymmetric Jacobian actions include thermal-softening, viscous-heat
and conductivity cross derivatives. Native cyclic positive spring solves
precondition the frozen momentum and conduction blocks; there is no joint-SPD
claim or new dense solver. Bounded native Newton/FGMRES, a whole-call budget
on joint source-state preparations and independent actual-state velocity and
cell-joule residual checks gate publication. Seed and transport work retain
their separate native bounds. Krylov storage is O(N*restart), and the original
source/preconditioner arrays are reused.

Unsupported algebraic source trials may backtrack within original support.
Accepted-source failures, exhausted work, arithmetic or energy refusal and
cancellation roll back the whole call. The native geometry-stop convention
still retains only completed coupled substeps. Explicit viscosity and
conduction CFL bounds are omitted from the stiff stage; transport, capillary,
geometry, source-support and caller maximum-step limits remain active.

Physical endpoint viscous heat enters internal energy. Backward-Euler kinetic
damping is accounted separately as `0.5*rho*pi*dz*sum(M_face*(du)^2)` and never
deposited as heat. The report identifies the actual first-order method and
retains thermal, mechanical and combined signed residuals. Stable large
implicit steps can underresolve fast physical heating; Newton convergence is
not a time-accuracy certificate, and timestep refinement remains necessary.

Four focused regression groups cover coupled first-order refinement, stiff
original rheology/conduction beyond the explicit diffusion bound, analytic
cross-block derivatives, source-trial backtracking and atomic late refusals.
Independently executed coupled Jacobian checks gave relative errors below
2.4e-10; timestep-halving errors fell by factors 0.5008--0.5085. A separate
stiff BE reference took a step about 771 times its initial explicit combined
diffusion limit, with 0.00342579827 J physical heat and 0.19290853318 J numerical
kinetic damping. These independent equation checks are distinct from native
Rust execution, which remains blocked by the workspace dependencies/storage.
All nine changed Rust sources pass pinned rustfmt parsing.

### Heat-aware adaptive implicit jet evolution

`ThermalCarreauJet::advance_implicit_adaptive`, under
`implicit-enthalpy-coupling`, compares exactly one original joint IMEX-Euler
step with two original joint half steps. The fixed and adaptive entries share
the same single-step kernel and final publication path. Only the two-half-step
trajectory is retained, without Richardson extrapolation or a new order claim.
`thermal::implicit::adaptive` exposes physical tolerances and the report.

Local error indicators cover area, radius, actual face velocity, original
caloric temperature and specific energy, per-cell physical viscous heat, and
the noncancelling sum of absolute cell-heat discrepancies. Relative scales use
changes over the trial interval, so a large energy datum, absolute temperature
or bulk velocity cannot loosen the local-change allowance. Heat scales use
only that interval's physical production.

A separate accepted-prefix account enforces `D <= D_abs + r*(H+D)`, with
`0 <= r < 1`, physical viscous heat H and numerical BE kinetic damping D.
Direct trial costs and double-double accumulated totals prevent a small
positive loss from disappearing beside an exhausted large history. Exactly
zero new loss remains admissible and allows growth during a quiet tail.
Absolute-only budgets also reserve future credit using current centered
kinetic plus nonnegative surface energy, releasing that scheduling reserve
when centered kinetic energy fits the remaining allowance. This heuristic
does not replace the exact prefix guard or claim a future capillary bound.

All coarse and fine trials, including discarded attempts, consume the shared
substep, original Newton and source-preparation budgets. Coarse/rejected heat
never enters the accepted cell or aggregate accounts. Nine reusable O(N)
checkpoint arrays supplement the original bounded Krylov workspace; no
trajectory archive or per-trial clone of the simulation is added. Original
transport/capillary CFL, geometry and material support remain authoritative.
Source, work, refinement, representation or conservation refusal and
cancellation preserve the complete pre-call state. A native geometry stop
publishes only earlier completed comparison intervals.

These are local solution/heat indicators and an actual interval-wide numerical
damping allowance, not a certified global physical-heat error bound. The
default entry retains first-order integration. Numerical damping is never deposited
as physical heat; explicit transport and nonlinear residual errors remain.

Focused regressions exercise stiff heat recovery and budget refinement,
accepted-only spatial heat, original source retention, late refusal/retry,
cancellation and exhausted-budget sub-ULP costs. Independent coupled
Carreau/WLF, nonlinear-Cv and conduction calculations over 0.02 s give
0.196349540849355 J physical heat in the resolved continuous reference.
Absolute damping allowances 0.002 and 0.0005 J give physical heat
0.194350019873239 and 0.195849542010369 J, with numerical losses
0.001999520976158 and 0.000499998839052 J. Heat, temperature and deformation
errors improve by factors about 0.25 when the allowance is quartered.
These independent numerical calculations and pinned rustfmt parsing do not
replace native Rust compilation/tests, which remain unexecuted in the current
dependency/storage-constrained workspace.

### Second-order energy-consistent thermal jet evolution

The method-aware entries `advance_implicit_with_method` and
`advance_implicit_adaptive_with_method` additionally admit
`ThermalJetImplicitMethod::StrangMidpoint`. Existing entry points retain
`ImexEuler` behavior. The second-order method composes an SSPRK2 conservative
transport/capillary half step, one original jointly solved implicit midpoint
viscous/conductive step, and a second SSPRK2 transport half step. Order two
applies to smooth admitted semidiscrete solutions with resolved nonlinear
error; limiter changes and material/source boundaries are not smoothness
guarantees.

The midpoint equation is the existing joint BE equation at half duration.
It solves conserved momentum and extensive caloric energy together, with
the same original rate/temperature viscosity, conductivity, analytic Jacobian,
bounded Newton/FGMRES and native block preconditioner. The actual midpoint
temperature is the original inverse of midpoint energy, not an arithmetic
temperature average. Endpoint momentum and energy are reflected about that
midpoint using differences from the transported source state. Both the
midpoint and reflected endpoint must remain inside original source support.
A supported midpoint never authorizes clipping an unsupported endpoint.

At the fixed areas of the dissipative subflow, the midpoint kinetic work
identity is `Delta K = h*<u_mid,Gq_mid>`; the same original cell power
`h*P_mid` enters internal energy and the cell/aggregate physical heat
accounts. Reciprocal Fourier exchange cancels across internal faces.
This subflow has zero BE kinetic damping. Actual reflected-state momentum
and cell-energy residuals are checked against the caller's physical
tolerances, after requesting half those allowances from the midpoint solve.
Explicit transport/capillary and nonlinear/roundoff errors remain in the
signed mechanical and combined reports. Zero implicit damping is not a
certificate of full-flow energy accuracy or passivity.

Midpoint is A-stable but not L-stable: an unresolved rapidly decaying mode
can reverse sign, and a large conduction step can exchange hot/cold states
while conserving total energy. The method-aware adaptive entry therefore
retains original state, temperature, specific-energy, cell-heat and
noncancelling total-heat comparisons. Its fine-path local estimates divide
coarse/fine differences by three and use cube-root scheduling. Only the two
fine steps are accepted, without Richardson extrapolation. The BE numerical
damping guard remains available for BE; a zero midpoint damping account
cannot replace the trajectory/heat tolerances.

One reusable O(N) checkpoint of area, momentum and energy restores the
beginning of a split if its later transport phase refuses a CFL or geometry
trial. Original source, work, residual, cancellation and whole-call rollback
semantics remain. Physical midpoint power is captured before reflected or
transport endpoint queries overwrite source scratch. Coarse/rejected heat
is discarded; the two accepted adaptive halves share the native reset/add
heat buffer. No second constitutive solver or per-trial simulation clone is
introduced.

Five focused regression groups cover second-order original nonlinear-caloric
refinement with and without conduction, a supported midpoint whose reflected
temperature leaves source support, late work/cancellation and clean retry,
geometry stopping after an incomplete split, and adaptive stiff-transient
recovery. Independent coupled equations at h=0.001 s give 0.103642992530559 J
physical heat with conduction, versus the resolved continuous value
0.103641958008171 J; all five state/heat errors decrease by about fourfold on
timestep halving. A stiff single 0.02 s midpoint step conserves energy while
producing only 0.01359096426 J and reversing velocity. The adaptive independent
calculation recovers 0.196349540849376 J against 0.196349540849355 J continuously,
with 21 accepted and 9 rejected comparisons. At an earlier transient,
tightening active local tolerances tenfold improves deformation, velocity,
both temperatures and physical heat by factors about 0.23--0.25. These are
executed independent equation checks. Pinned rustfmt parses the changed Rust
sources; native compilation and Rust regressions remain unexecuted under
the current dependency/storage constraints.

The construction uses the second-order half-step composition of Strang
(SIAM J. Numer. Anal. 5, 1968, DOI 10.1137/0705041) and the quadratic-work
identity of implicit midpoint (Hairer, Numerical Geometric Integration,
1999, section III.3). Their order and algebraic identities do not establish
unconditional source admissibility or resolved stiff dynamics.

## Adaptive source-owned turbulent gas heat and species transport

`FvFlow::advance_turbulent_gas_species` advances a complete physical window
through the original material network's coarse-versus-two-half-step controller.
Native geometry, accepted carrier masses and momentum `nu_t/Sc_t` face
conductances remain fixed for that flow epoch. Every probe refreshes actual
composition, donor partial enthalpy and the original caloric/EOS/transport
endpoint. Only accepted fine steps contribute to the compensated stream ledger.
The result `TurbulentGasSpeciesAdvance` retains local species/temperature error
ratios, bounded work and accepted-step statistics, original-to-final physical
audits and carrier drift against the same original accepted masses.

`diffuse_turbulent_gas_transport` and `advance_turbulent_gas_transport` add
turbulent heat to the snapshot and adaptive entries respectively. Source-derived
`k_t = rho cp(T,Y) nu_t,heat / Pr_t` uses the retained thermal eddy field and
positive Prandtl number; momentum `nu_t/Sc_t` remains the independent species
law. Graded half-cell series resistances give antisymmetric internal heat
`Q = dt K (T_left - T_right)`. Each probe samples heat and species from the
same original state and sends both to one material reconstruction. Donor partial
enthalpy and formation references are carried once; additional heat has its own
accepted-only signed ledger. Solid, exterior and zero-conductivity faces supply
no heat in this component. Molecular conduction and external wall heat retain
their separate owners.

The original 256-parcel and 4096-stream bounds apply; heat and species share the
stream/whole-window report bound. Typed withdrawal and caloric-bracket overshoot
may shorten adaptive trials. Other source-domain errors are fatal. Local-error,
step and evaluation limits bound refinement; rejection/cancellation never
publishes a partial candidate. Fixed whole-window mass, each-species, actual
source-relative energy and endpoint entropy budgets are never multiplied by
the number of accepted steps. The source-relative energy audit removes common
formation references before arithmetic. Carrier drift is also measured against
the initial flow state, not reset after every substep.

Snapshots and the retained adaptive fine path remain first order in time.
Local control is not a global-error or turbulence-model validation certificate.
The held eddy/geometry epoch does not refresh pressure or momentum; subsequent
flow consumers must use the returned material endpoints. The original planar
or finite-volume turbulence closure limitations remain.

Focused regressions are authored for analytic species relaxation, refreshed
partial enthalpy, typed withdrawal recovery, variable-heat-capacity heat
exchange, formation-reference invariance, an independent combined RK4
reference, accepted-only heat accounting, independent Prandtl/Schmidt response,
thermal-bracket recovery and atomic late work/source/cancellation failure.
Scoped pinned formatting and whitespace checks passed. Full native integration
execution remains blocked by the pinned asupersync dependency exceeding this
environment's build resources before the changed crates compile.


### Partial enthalpy on the same constituent diffusion transfers

`MaterialAdvectionStep::diffuse_species_with_enthalpy` adds thermal storage to
the original common-diffusivity owner. One callback supplies the original
partial specific enthalpies in J/kg at each conducting face, frozen for the
interval. Every actual signed constituent transfer contributes `sum(h_k dm_k)`;
one cell loses exactly the energy its neighbor receives. The existing positive
mass relaxations, graded-grid conductances, ordered sweeps and mass/partition
budgets remain the only species solve. Source face queries occur once before
the sweeps, and both inventories publish together only after final admission.

Because the barycentric diffusive carrier flux is zero, a shared enthalpy datum
is subtracted before products. Species-specific formation/latent differences
remain. Error-free subtraction and product transforms retain their low parts
until opposing species contributions cancel; compensating already-rounded
products alone is insufficient. Compensated cell increments update `m*h` storage
in centered coordinates; an absolute joule budget includes local unrepresented
work, closed-domain storage/transfer defects, and rounded carrier-flux defects
weighted by the partial-enthalpy range. The endpoint audit also retains low
parts of mass times the represented enthalpy change; opposing local storage
errors cannot disappear through rounded products or a zero global balance. Large formation-energy totals cannot
enlarge that budget. Thermal face/scratch/output storage is included in the
original allocation ceiling. Source error, cancellation or failed mass/heat
admission returns no partial candidate.

The returned endpoint enthalpy already contains the reported diffusion energy.
Consumers must invert it with the returned composition using their original
caloric law, re-evaluate face properties as needed, and accept both inventories
together. This component supplies no EOS, source identity, temperature solve,
reaction or pressure-work source, arbitrary multicomponent molecular law, Soret
term or coupled temporal-order guarantee. Freezing partial enthalpies is a
declared approximation for a changing state. Focused G0/G1/G3 regression
definitions cover analytic binary relaxation/energy flux, multicomponent
isothermal caloric inventory, common-datum invariance, retained small heat amid
opposing formation terms, unchanged mass-only results, and unrepresented-work,
source, storage and cancellation refusals. These
new regressions require native execution; the authoring environment went
offline before that check completed.

## Original gas calorics on constituent diffusion

`GasSpeciesEnthalpySource` binds the full canonical species bank, original
NASA-9 formation references, pressure context and transport data to the native
species/enthalpy diffusion owner. Every cell preparation converts the actual
constituent masses to molar amounts, checks their partition of the immutable
carrier to machine roundoff, and inverts the original composition-dependent
enthalpy law. The same resolved state supplies temperature, Cp, conductivity,
viscosity and EOS density. Zero local abundance retains the full source bank.
The full bank must admit a continuous increasing caloric branch, including
species absent from a particular cell; discontinuous seams do not get hidden
by a diffusion step.

Face partial specific enthalpies come from each original species law at the
physical face temperature, interpolated with the actual graded cell widths.
`GasSpeciesDiffusivity::Binary` reuses an explicit `BinaryDiffusionMaterial`
for exactly its two species, refreshing its m²/s law at the reconstructed T/p
while preserving auxiliary source coordinates and property selection. The
`Frozen` variant retains an explicitly authored common mixing field, including
zero-diffusion barriers. Viscosity never substitutes for molecular diffusivity,
and one binary coefficient cannot silently represent a larger molecular mixture.

`diffuse_thermochemical_species` advances the existing positive species kernel
and its exact constituent-carried energy account, then reconstructs the endpoint
from both returned inventories. Initial/endpoint caloric residuals and the
original transport defect share an absolute joule budget. Owned prepared arrays
are reserved from the native storage ceiling before allocation; immutable source
banks and temporary original material-evaluation scratch remain additional.
Source, identity, domain, storage or cancellation refusal returns no partial
candidate. The returned h already contains diffusion energy.

This entry freezes D and face partial h over its declared interval; repeated
steps must refresh them. It is an isobaric fixed-carrier constitutive closure:
EOS density is returned for subsequent flow/pressure coupling and never replaces
the accepted mass. It introduces no pressure work, reactions, arbitrary
multicomponent molecular friction, Soret term or coupled temporal-order claim.
Regression definitions cover variable-Cp analytic binary exchange, formation
reference invariance, graded physical face interpolation, absent species,
source-dependent diffusivity, validity axes and atomic refusals.

## Enthalpy through original species reservoirs

`diffuse_species_with_enthalpy_and_boundaries` observes each actual signed
constituent transfer from the existing positive reservoir/interior maps.
One original partial-specific-enthalpy vector in J/kg is supplied per active
boundary face. The same local face law applies to counterdiffusing species;
the common carrier flux is zero. The original film/half-cell conductance,
ordered sweeps, positive species inventory and mass/partition gates remain
authoritative. Empty boundaries reproduce the closed enthalpy entry exactly.

Interior and boundary energy share the original centered partial differences
and error-free products. A boundary credits the same retained energy components
to its cell and the independent external account once. The shared endpoint
audit compares represented storage and transported energy with external heat,
rather than incorrectly requiring a closed-domain zero. Per-cell storage
mismatches and every boundary mass-frame defect still cost absolute joules.
Internal cumulative energy and species receipts retain low parts; public
scalar readbacks round only at the reporting boundary.

Per-face and aggregate receipts retain signed energy and positive/negative
NET map energy after cancellation among species. Large opposing reservoirs
or a reversal during an interval cannot disappear behind net domain heat.
This energy is already contained in the returned enthalpy. The generic
reservoir observer forwards both kinds of transfers; the mass-only consumer
is a zero-sized no-op. Active thermal face caches, retained accounts and
readbacks are admitted against storage before allocation, and source refusal
or cancellation publishes no candidate.

`GasSpeciesEnthalpyBoundarySource` binds the exact canonical reservoir list to
explicit PHYSICAL face temperatures and the original gas source. Its partial
enthalpies use the same NASA evaluations as interior faces, including species
absent from a cell. `diffuse_thermochemical_species_with_boundaries` resolves
the endpoint composition and total enthalpy together, refreshes gas transport,
and charges initial/endpoint caloric residuals alongside native open-energy
defects. Species-specific formation shifts transform external chemical energy
and stored enthalpy together; a common mass-specific datum cancels from the
zero-carrier diffusion flux.

Reservoir composition and a finite mass-film resistance do not prescribe a
thermal interface temperature. The supplied face trace must come from its
own thermal boundary/interface owner; no donor-temperature, Fourier-film,
reaction, extra mixing-heat or pressure-work source is inferred here. The
reservoir is fixed and owns no evolved storage. Source coefficients are frozen
over this entry's interval; coupled temporal refresh has its own driver.

Focused regressions include exact finite-film exchange, simultaneous interior
and boundary transport, empty-boundary equivalence, signed/gross heat reversals,
common-reference invariance, retained formation and species-product remainders,
and source/storage/cancellation refusals. Original gas integration fixtures
exercise an explicit face temperature distinct from cell temperature and
formation-reference covariance.

## Adaptive species and enthalpy with original material feedback

`SpeciesEnthalpyDiffusionModel` provides one original state-dependent caloric
and transport closure. `diffuse_species_with_enthalpy_adaptive` prepares the
original state, advances a positive half-step predictor, and prepares its
midpoint properties. A full-step corrector then starts from the original
substep inventories using those properties. Predictor heat is never deposited
again. The endpoint source image supplies the next substep and is returned
as `state` on acceptance; a gas caller obtains refreshed Cp, conductivity,
viscosity, EOS density and diffusivity without repeating the final inversion.
`GasSpeciesEnthalpySource` implements this same generic closure directly.

Complete intervals at n and 2n physical refresh steps estimate smooth
second-order time error through |fine-coarse|/3. Temperature in K, specific
enthalpy in absolute J/kg, and each constituent mass fraction have independent
allowances. Formation-energy totals cannot enlarge a thermal tolerance.
The accepted candidate is the original positive fine trajectory, with no
extrapolation, clipping or composition normalization. The prior fine solution
becomes the next coarse comparison. This is estimated temporal accuracy,
not a certified global error bound; sharp constitutive changes and phase corners
need their own source validity and resolution controls.


Complete relaxation can make coarse and fine endpoints agree after both have
skipped a thermal transient. The driver therefore compares sampled D and
error-free centered species partial-enthalpy differences at the start and
midpoint. For changing coefficients, the accepted outer face relaxation must
not exceed one half. That rate is recovered over the full physical refresh
step, so increasing the inner frozen sweeps cannot bypass the check. Exact
unchanged coefficients retain the large-step path. This anti-saturation check
uses counted face queries without another source inversion or face buffer;
it supplements the temporal estimate rather than certifying a constitutive law.

Every predictor, corrector and discarded trajectory contributes to total
pair/component work and source-query counts. The storage gate reserves live
comparison states and three source preparations before allocation; obsolete
coarse source images are dropped before finer work. Original source-table
storage and temporary source-evaluation scratch remain additional. A refused
source state, depleted budget or cancellation publishes no candidate.

Accepted physical correctors accumulate their original constituent-mass and
enthalpy defects. A separate initial-to-final audit uses error-free subtraction
and products for represented storage. Internally, the frozen kernel returns its
existing compensated energy increments by move; both components survive
substep composition, even when a public scalar heat readback rounds away a
small term beside formation energy. Initial/final caloric-inversion residuals
remain explicit in the joule budget, and the largest accepted-stage residual
is reported separately. Rejected trajectories consume resources but are not
added to the accepted physical history.

Focused original-kernel regressions cover the isothermal frozen limit,
independent variable-Cp analytic convergence, a smooth nonlinear heat-of-mixing
case that defeats unguarded step doubling, simultaneous thermal/diffusivity
feedback against an independently refined ODE, source/work/storage/cancellation
refusals, and exact conservation of small heat between large opposing formation
references. A separate integration regression couples the original NASA source,
binary material diffusivity and the native carrier to an independent two-cell
mass/energy ODE. These remain fixed-carrier constitutive diffusion steps;
pressure/volume evolution, reactions and nonideal molecular cross-effects
require their original physical owners.

## Local SGS with general prescribed moving walls

`Grid3::prepare_sgs_with_moving_walls` and the corresponding serial/pooled
step entries use local Smagorinsky, WALE or Vreman closures on selected fixed
halfway-wall links, including interior solids and spatially varying tangential
drives. The complete wall-link field remains with the native moving bounce-back
owner. SGS and kinematics share each missing Cartesian neighbor's actual axial
wall-link velocity; diagonal-only contacts do not invent an axial wall.

The ordered link mapping uses fixed storage and no per-cell probe search or
additional field allocation. Existing collision-family rates, supplied molecular
viscosities, Mach/admission checks and eddy-to-molecular limits remain intact.
The collision field is prepared from fresh molecular values on each call.

`ForcedGrid3` exposes the same preparation, held-force steps and endpoint-force
advances. The accepted physical velocity uses the current `F_n/2` correction
before gradient fitting; `F_(n+1)` admission and publication remain in the original
force owner. Native population, force, clock and wall-load transactions are reused.
SGS diagnostics report instantaneous modeled shear transfer and add no heat.

Regression groups are provided for analytic wall gradients on both lattices and
all three models, zero-drive/disabled controls, serial/pooled/native equivalence,
force histories, cancellation and refusal. This remains the exploratory local
closure on fixed staircase geometry; it does not add wall functions, geometry
remapping, a turbulent-flow qualification or an integrated energy balance.


### Molecular stress with phase resistance and changing density

`advance_viscous_phase_resistance_thermal_volume_momentum` composes the
existing full-stress molecular owner with source-resolved mushy resistance,
physical heating, thermal volume flow and conservative M/Mh/P transport.
Original material receipts supply molecular viscosity; actual M/V supplies
inertia. Darcy resistance resolves the accepted endpoint phase on the same
material chart. Molecular and Darcy heat each enter caloric storage once;
numerical losses, reservoir work and source-pressure work retain their separate
accounts. Noncancelling caloric and mechanical defects share the existing
energy budgets. Source, solve, cancellation or budget refusal publishes no
partial inventory.

This is a first-order split with frozen original molecular viscosity and
endpoint Darcy coupling. It provides neither simultaneous nonlinear shear
coupling nor exact solid locking, moving skeletons or free surfaces. Initial,
incoming, heated and transported phases must remain mobile.

Four new focused regressions cover closed kinetic/caloric balance with both
heat sources, unequal phase densities and volume response on every axis, exact
fully-liquid parity with the original viscous route, and atomic source/solid/
cancellation refusals with retry. These and three existing phase-resistance
tests passed in a native production-source harness; this is not a full workspace
build or a temporal-convergence certificate.

## Original material rheology with curved-wall SGS

`MaterialScaling::step_curved_walls_with_material_rheology_and_sgs` couples
the original source response to the native off-lattice wall owner. Each fluid
cell is sampled once at incoming physical deviatoric shear. One accepted
Guo-centered full gradient supplies both molecular rheology and the explicit
Smagorinsky, WALE or Vreman closure. Selected axial wall links use their actual
intersection distances and supplied wall velocities; diagonal intersections
retain the original BFL population/load geometry.

The exact same frozen wall velocities drive gradient, population exchange,
impulse, torque and wall work. Original local molecular viscosity supplies the
eddy/molecular ratio gate and molecular diagnostic. Collision family, TRT
product, independent bulk and higher-order controls remain unchanged. There
is no second gradient, full-grid clone, implicit material solve or heat deposit.
The original resolved-shear nonlinear-material approximation and fixed
reference-density/scaling limits remain explicit.

This reuses `BoundaryRheologySgsWorkspace3` and its two source buffers.
Applied material samples and the SGS/selected-wall report publish only after
the original curved collision/stream/load transaction succeeds. The complete
returned native receipt also retains other walls, body force, BFL mass change,
open and periodic transport. Refusal/cancellation leaves accepted populations
and material samples unchanged.

Three focused regressions are authored for affine/quadratic wall profiles on
every axis and all five collision families, independent material/eddy
viscosities and momentum/angular-momentum/work accounting, exact constant-source
native recovery, single source reads, and late-refusal replay. Scoped pinned
formatting and whitespace checks passed; native execution remains pending the
pinned dependency build's resource requirement. This composition adds no
turbulent-flow, unresolved roughness or wall-resolution validation claim.

## Source molecular conductivity in turbulent gas transport

`FvFlow::diffuse_turbulent_gas_transport_with_molecular_heat` and
`advance_turbulent_gas_transport_with_molecular_heat` explicitly compose
original gas molecular conduction with the existing turbulent heat/species
update. Each current material state supplies
`k_total = k_source + rho cp nu_t,heat / Pr_t`. These total cell conductivities
enter the graded half-cell series resistance; independently harmonic molecular
and eddy face conductances are not added. Molecular heat therefore persists
when either cell's thermal eddy viscosity vanishes.

Every adaptive probe refreshes original conductivity, density, heat capacity
and donor partial enthalpy at its own state. One original caloric update receives
all heat and species streams. Its accepted-only heat ledger reports total
molecular-plus-turbulent paired heat; callers must not add another molecular
conduction step for the same interval. Existing turbulent-only entries retain
their original component semantics. This supplies neither molecular species
diffusion nor Soret/Dufour, external-wall or reservoir transport.

The original bounded first-order snapshot/adaptive owner, formation-reference
invariance, fixed mass/energy/entropy allowances and cancellation behavior
remain. Focused regressions cover the zero-SGS molecular limit, total half-cell
resistance, source retention and adaptive accuracy. New numerical code/tests
are scoped-formatted and statically reviewed; native execution remains blocked
before the changed crate by the pinned dependency build's resource requirement.

## Transient native thermal boundary inputs

The existing `ThermalOpenField3::set_reservoir_temperature` supports changing inlet and backflow reservoir temperatures without reconstructing the scalar field. `set_wall_boundary(face, law)` now updates an existing external stateless wall uniformly over its original face; prescribed temperature, signed inward heat flux, ambient-film transfer/data and insulation can change between accepted intervals. `set_internal_wall_boundary(law)` changes only the fallback law on unassigned fluid-to-solid links. Explicit immersed bodies retain their link ownership, stored heat and conduction graph. These setters preserve native face geometry, scalar populations, flux history and accepted time; they cannot create or replace finite thermal storage.

`set_finite_wall_heat(face, heat_per_step)` changes the external source of each retained finite wall cell in original x-fast surface order. The field covers the complete wall array, is finite, and is zero in slots adjacent to solid grid cells. Negative source values represent cooling. Stored temperatures, capacities, contact resistances and prepared wall-conduction factors remain unchanged; no energy is deposited by the setter. The original finite-wall solve credits each active source once and balances it against fluid heat, solid storage and reciprocal wall conduction. `wall_states` exposes current sources, while `boundaries` retains the original finite-wall declaration.

All inputs are admitted before mutation. A refused update leaves both inputs and physical state unchanged. A subsequent failed or cancelled step retains the newly authored inputs and the prior physical epoch, matching reservoir-update semantics; the next attempt rebuilds scratch from the retained state. The Boussinesq owner forwards all three APIs, pre-admits prescribed wall temperature through its density-departure envelope, and retains the original actual-trace checks for flux and ambient-film laws. Native held-force, varying-force, scalar-diffusion and tensor-diffusion steps use these same boundary owners without a second heat solve or extra source allocation.

## Spatial-force and buoyant open-boundary SGS

`ForcedOpenBoundaryGrid3` now prepares and advances the existing local
Smagorinsky, WALE and Vreman closures with each cell's accepted spatial Guo
half-force. Held-force and next-force entries share the native reservoir/wall
gradient and original collision/stream transaction. SGS samples only the
opening physical velocity; the supplied endpoint force retains its separate
native reconstruction and trapezoidal impulse semantics.

`BoussinesqOpenGrid3::step_serial_with_sgs` and `advance_serial_with_sgs`
compose that same momentum closure with `alpha_t = nu_t / Pr_t` and the
original D3Q7 molecular thermal relaxation. Fresh molecular collision/thermal
fields are explicit, augmentation is frozen for the interval, zero SGS retains
their exact bits, and eddy/molecular limits refuse rather than clip. One shared
coefficient preparation serves passive and buoyant thermal owners without
additional per-step arrays or population clones.

The original joint owner retains heat-source half-kicks, reservoir/backflow
exchange, finite conducting walls, temperature-dependent endpoint buoyancy,
and separate body/background/buoyancy impulse accounts. Successful publication
alone advances flow, heat, walls, force epochs and clocks or exchanges a supplied
endpoint-background buffer. Late source, ratio, force or cancellation refusal
does not publish a partial physical state. Modeled shear transfer remains
diagnostic; neither it nor mechanical-force work is added to heat.

Four focused regression groups cover analytic spatial-force gradients and
native composition, exact zero-SGS trajectories, active buoyant feedback with
independent fluid/solid heat and impulse accounts, and warmed late-refusal
replay. The first three exercise D3Q19 TRT and D3Q27 cumulant lattices with all
three closures. Scoped pinned formatting and whitespace checks passed before the execution
environment went offline; the final source was recovered from those reviewed
edits through GitHub. Native execution remains pending the pinned dependency
build's resource requirement.
This is the existing explicit fixed-reference-capacity Boussinesq split with
local exploratory closures, not a variable-density total-energy model,
buoyancy-production turbulence closure or new temporal-order qualification.

## Forced open-flow SGS momentum and heat

`ForcedOpenBoundaryGrid3::step_serial_with_thermal_sgs` and
`advance_serial_with_thermal_sgs` combine explicit local Smagorinsky,
WALE or Vreman momentum/heat transport with held or advancing spatial force
fields on D3Q19 and D3Q27. Preparation uses the accepted physical velocity,
including each cell's Guo half-force correction F_n/2; a proposed endpoint
force does not enter that incoming fit. The existing shared coefficient
preparation supplies one eddy-viscosity field to native shear relaxation
and to D3Q7 heat diffusivity through `alpha_t=nu_t/Pr_t`.
Callers supply fresh molecular transport each time, preserving the original
closure, coefficient-ratio and disabled-model admissions.

The original joint forced-thermal transaction retains force-buffer exchange,
reservoir transport, prescribed thermal-source half-kicks, finite stationary
wall storage, wall-work receipts, histories and counters. One existing SGS
workspace supplies all coefficient storage; the wrapper introduces no
per-step allocation or alternate physical state. Publication occurs only
after the original native transaction succeeds. Refusal or cancellation
leaves accepted physical fields and caller-owned endpoint forces unchanged;
scratch may change.

This serial API remains fixed-reference-capacity temperature transport:
modeled mechanical transfer and body work are not deposited as heat. It
does not add a total-energy model, wall function, pooled force/thermal owner
or empirical turbulence validation.

Three regression groups are authored for both lattices and all three
closures: constant physical velocity under spatial forcing versus a
wrong-force-frame negative control; coefficient/Prandtl agreement, native
joint-transaction equivalence, disabled controls and storage reuse; and
invalid endpoint force, late thermal failure, diffusivity-ratio refusal,
cancellation and retry. Pinned-nightly syntax/format checks and source
review passed. Native compilation and execution of this extension remain
unverified: real-graph builds exhausted memory upstream in `asupersync`
before reaching LBM. No test-pass or performance claim is made for those
unexecuted cases.

## Analytic full-gradient local SGS response

`sgs3::eddy_viscosity_tangent` returns the original Smagorinsky, WALE or
Vreman viscosity plus all nine analytic derivatives with respect to
`G[i][j]=d(u_i)/d(x_j)`. `vreman_eddy_viscosity_tangent` additionally
retains three directional filter widths, assigned to derivative coordinates.
Both return `SgsViscosityTangent { viscosity, derivative }`.
The original primal closure supplies the viscosity bit-for-bit. Fixed-size
invariant differentiation introduces no allocation or production finite
differences. WALE retains the full gradient's strain and rotation dependence;
Vreman retains compensated weighted minors and exponent-scaled arithmetic.

Derivatives have dimensions of length squared, with coefficients and filter
widths held fixed. A stress Jacobian must include BOTH terms of
`d(2 nu_t S)=2 nu_t dS+2 S(dnu_t/dG:dG)`; the returned viscosity derivative
alone is not the stress Jacobian. Overflow or underflow of an otherwise
nonzero derivative explicitly refuses instead of clipping.

At Smagorinsky zero strain, all-zero WALE gradient and Vreman zero/rank-one
states, the response uses a documented zero generalized-derivative
selection. These are not claims of classical differentiability: a retained
rank-one Vreman example has nonzero one-sided derivatives. At nonzero WALE
gradients with zero traceless symmetric square, the ordinary derivative is
zero. No new nonlinear convergence, turbulent-flow validation, wall-function
or transition claim is introduced.

Twenty native numerical tests passed on pinned nightly: ten original
closure tests plus ten tangent groups covering all nine entries against
independent centered differences, exact primal equality, isotropic rotation
covariance, parity and homogeneity, unit/filter scaling, WALE near-wall cubic
behavior, anisotropic Vreman axes/aspect ratios/rank limits, complete modeled
stress response, disabled models and unrepresentable derivatives. The build
compiled the complete original `fs_math` crate and the actual
`model.rs`, `vreman.rs` and `tangent.rs` sources without runtime
adapters, numerical stubs or substituted kernels. This establishes the
numerical slice, not full `fs-lbm` or downstream solver integration.

## Implicit conjugate heat with molecular and phase resistance

`FlowField::advance_isotropic_implicit_conjugate_viscous_phase_resistance_thermal_volume_momentum_with_baths`
combines the native full molecular stress stage, implicit spatial conduction,
finite solid storage, bath/contact exchange and phase-dependent Darcy drag on
the original material chart. Molecular inertia uses actual M/V; original checked
phase sources supply its viscosity once. Its physical heat joins authored heat
once, then the existing resistance owner iterates Darcy heat, implicit thermal
state and density-compatible volume flow from the same original fluid and solid
inventories. Only the accepted trial transports M, M*h and P and publishes a
finite-wall endpoint. Fixed-solid mass and momentum stay outside fluid routes.

Source conductivity uses the native analytic enthalpy Newton/Krylov graph.
Whole-call residual-evaluation and Krylov ceilings include rejected outer
trials. Molecular assembly, implicit heat, conservative caloric transport and
Darcy heat/constitutive defects enter the shared caloric and both mechanical
budgets exactly once. Physical molecular/Darcy heat, bath/contact transfers,
pressure/source work and numerical mechanical losses remain separate accounts.
Cancellation or any source, heat, mechanical or work-budget refusal publishes
no partial momentum, material inventory or heated wall.

G1/G3 checks independently reconstruct finite-wall/contact Fourier heat,
unequal-density latent volume expansion and combined fluid/wall/advection
mass/enthalpy balance while both physical dissipations are active. Zero Fourier
exchange preserves the prior composed molecular/resistance endpoint exactly.
G4 checks cancel after a completed heat trial, exhaust aggregate residual work,
and reproduce the unchanged accepted result on scratch reuse. Focused execution
uses actual material, native sparse, deterministic reduction, LU, FGMRES,
enthalpy/conduction and flow sources; full Cargo/execution-context integration
remains a separate verification gate.

The split retains original-state molecular viscosity, implicit pre-advection
conductivity and transported-endpoint phase resistance. It is first order and
requires refinement; it does not infer simultaneous endpoint shear feedback.
Fluid phases must remain mobile. Exact solid locking, moving skeletons and
free surfaces remain outside this owner.

## Thermal viscosity feedback inside evolving mushy flow

`FlowField::advance_isotropic_implicit_conjugate_thermoviscous_phase_resistance_thermal_volume_momentum_with_baths`
returns `ImplicitThermoviscousResistiveThermalVolumeMomentumStep` and closes thermal
shear and independent bulk viscosity inside the original combined molecular,
implicit heat, phase-density and Darcy timeline. Each Darcy heat trial reuses
the native thermal-viscosity fixed point on ORIGINAL fluid M/P/h and held-solid
inventories. Exact native pretransport caloric samples resolve molecular source
receipts; original M/V supplies inertia even when the heated chart density has
changed. In particular, converting a source kinematic viscosity with the heated
source density must not replace original mass in the molecular kick.

The existing mobile phase and source-resistance admission checks remain mandatory:
chart and material identity, exact phase, exact source density, and the original
mushy volume/permeability evidence must all agree. The native coefficient-only
stress seam consumes the admitted shear and independent bulk response with
`2 mu dev(D) + zeta div(u) I`. The converged original-inertia molecular momentum
then enters the SAME Darcy pressure, volume-source and conservative M/Mh/P owner.
Darcy resolves at the transported endpoint, while molecular viscosity and
conductivity resolve at the pretransport thermal endpoint. This is a first-order
mechanical split requiring refinement, not a simultaneous unsplit flow endpoint.
Rate-dependent constitutive feedback is outside this opt-in.

Both physical heat sources enter the same native implicit interval once. The
reported used molecular heat retains its non-cancelling difference from actual
physical molecular heat. Heat assembly, that mismatch, implicit heat, caloric
transport and Darcy heat/constitutive defects consume the shared caloric and
both mechanical energy gates once. Pressure/source and boundary work, molecular
and Darcy heat, and numerical time/scatter/transport losses stay separate.
Accepted source, heat, momentum and finite-wall receipts belong to one trial;
rejected trials never advance physical inventory.

`ThermalViscosityFeedbackControls::maximum_trials`, the heat residual and Krylov
limits, and molecular `max_sweeps` bound the WHOLE call, including all nested
rejected Darcy trials. The report counts every source callback, including Darcy
endpoint queries. The original frozen-coefficient path retains its cached
reconstruction. Liquid-only thermal and thermal/rate entries retain their own
phase-support rules through the shared preparation loop. Exact immobile phases,
moving solid skeletons and free surfaces still require their existing owners.

G1/G3 native tests exercise unequal phase densities, source-density-dependent
molecular viscosity, independent bulk response, finite-wall/contact conduction,
physical molecular and Darcy heat, latent volume source, and conservative
fluid/wall/advection mass and enthalpy balance. An independent native kick at the
accepted source reproduces the original-M result and discriminates against
substituting hot source density into inertia. G4 checks stale source evidence,
cancellation after an actual heat solve, exhausted whole-call heat, molecular
and trial budgets, and exact retry from unchanged finite-wall history. Existing
frozen combined, liquid thermal-feedback and joint thermal/rate regressions pass
with the same factored preparation. Focused execution uses original physical
kernels; full workspace Cargo/execution-context integration remains a separate
verification gate.

## Joint thermal and endpoint-rate response in evolving mushy flow

`FlowField::advance_isotropic_implicit_conjugate_thermorheological_phase_resistance_thermal_volume_momentum_with_baths`
returns `ImplicitThermorheologicalResistiveThermalVolumeMomentumStep`. It extends
the native mobile-mushy thermal/Darcy composition with a true nonlinear molecular
endpoint solve at each exact native caloric state. The molecular resolver receives
the supplied phase and physical RMS `sqrt(2 dev(D):dev(D))` in 1/s, and returns
original source evidence plus independent bulk response. The existing coupled
momentum rows refresh those coefficients after every sweep; acceptance checks the
true nonlinear row residual and physical stress heat at the accepted endpoint,
not a small coefficient change. That physical heat closes the same original
implicit fluid/finite-solid thermal interval.

Original mass, momentum and M/V face inertia stay fixed across every nested
trial. Temperature- or phase-dependent source density belongs to the source
viscosity conversion and exact phase evidence; it never replaces original mass.
Exact phase, chart/material identity, source density and mobile law/transport
admission are shared with the existing native owner. Source admission is now
separate from forming a Darcy coefficient: a molecular law is not rejected
because its UNUSED `mu * inverse_permeability` would overflow. Real Darcy
consumers and the existing shared-source thermal-only entry retain that product
and its original representability gate.

The new API requires a separate `resolve_resistance(cell, transported_phase)`
callback for the explicitly supplied pore-scale Darcy response. A macroscopic
molecular shear rate does not determine pore shear. The implementation never
fabricates a zero or frozen molecular rate for Darcy and never substitutes one
source family for the other. Both source families retain the same material chart,
current density and phase-support evidence; their accepted receipts remain
distinct. Molecular temperature/rate response closes before Darcy, pressure and
transport, while Darcy constitutive response closes at transported phase. This
remains the declared first-order mechanical split requiring refinement.

The generic native coupled owner carries the accepted rheological rate receipt
through its exact heat and Darcy trial without re-query. Whole-call source-field
updates, per-cell molecular plus Darcy queries, molecular sweeps, heat trials,
residual evaluations and Krylov work include discarded nested attempts. Physical
molecular and Darcy heat enter once; assembly, heat mismatch, implicit thermal,
transport and constitutive defects retain the original non-cancelling shared
caloric/work gates. Boundary work and numerical losses remain distinct. No new
SGS, moving skeleton, free-surface, exact-solid momentum or unsplit flow claim is
introduced. Existing liquid-only and thermal-only entries retain their scopes.

G1/G3 tests use a manufactured NONisothermal, density-changing mobile mixture,
a genuine molecular temperature/rate source surface and an independently
authored pore response. They check exact accepted source/rate agreement, true
native momentum residual, original-M nonlinear mechanics against a separate
native reference and a deliberately wrong hot-density inertia, physical heat,
finite-wall storage, compatible volume source and conservative fluid/advection
mass/enthalpy balance. A source-admission discriminator keeps a finite molecular
response despite an overflowing unused drag product, admits the independent
Darcy law, and verifies that real Darcy and legacy thermal-only product guards
still refuse. G4 checks late molecular and Darcy source refusal, cancellation
after actual heat work, whole-call work exhaustion and exact retry from unchanged
histories. Native focused execution also retains the existing liquid thermal/rate
and thermal-only mushy regressions; full workspace Cargo/execution-context
integration remains a separate verification gate.

## Actual-material-density SGS momentum and separated energy transfer

`FlowField::advance_material_momentum_with_viscosity_and_sgs` composes the
original source-viscosity predictor, force/pressure update and accepted
mass-route transport on closed fixed geometry. Original cell M/V supplies
every SGS dynamic-viscosity conversion and local molecular-ratio limit;
native half-cell masses still supply inertia. Source reference density does
not replace either authority. Smagorinsky, WALE and Vreman retain their shared
filter, native reciprocal-edge and explicit exploration policies.

One retained backward-Euler predictor refreshes SGS on its current native
velocity before the first solve and after every coupled sweep, before the
same-velocity residual. Original molecular viscosity is sampled once per
fluid cell and held for the interval. The added stress is deviatoric
`2 rho nu_t dev(D)`; the original molecular `2 mu D` response retains its
existing isotropic part. SGS supplies no bulk viscosity. Zero model coefficient
validates the configuration and then uses the original molecular path exactly.

Normal-cell work and each reciprocal shear-edge work use the SAME total
dynamic viscosity as momentum. At a material edge, total cell viscosities
enter harmonic transmission first. Resistance shares distribute the common
stress work, then each cell's original molecular and SGS fractions divide its
local loss. Wall traction and wall heat retain the existing molecular closure.
No independently averaged SGS stress, source-reference-density conversion or
second momentum solve is introduced.

`MaterialViscosityReport::cell_sgs_transfer_j` and `sgs_transfer_j` retain
modeled removal from resolved kinetic energy separately from molecular heat.
Only molecular heat enters the original specific-enthalpy inventory. The
existing kinetic/work audit includes both physical channels while retaining
backward-Euler loss, scatter loss and representation work separately. Source,
cap, solve, heat-insertion or cancellation refusal returns no partial candidate.

This remains a first-order split with frozen material density/viscosity during
the predictor. The modeled transfer is not an SGS energy-storage, cascade or
heating equation. Active pressure-mouth SGS remains outside this entry, as do
a variable-density total-energy qualification and resolved-turbulence accuracy
claims. Ten focused regression groups cover analytic graded-density decay,
deviatoric stress, harmonic work partition, molecular-only walls/heat, every
adjacent edge cap, zero-model identity and refusal/replay. Source parsing,
scoped formatting and independent arithmetic were checked; native execution
remains blocked before fs-lbm by the pinned asupersync dependency build.

## Native open-grid Boussinesq material transport

`BoussinesqOpenGrid3::material_transport_workspace` binds an original
`MaterialScaling` reference, physical dx/dt and an explicit signed Kelvin/theta
map to reusable collision, thermal-relaxation and query-temperature arrays.
The complete native scalar-temperature envelope must map to positive finite
Kelvin. Workspace and both returned source-sample vectors have explicit payload
budgets; original source owners' internal allocations are outside those caps.

`step_serial_with_material_transport` and
`advance_serial_with_material_transport` sample each original fluid material
once at opening and once at the actual candidate endpoint. The existing
`BoussinesqMaterial` admission bounds density/Cp departures and retains reference
rho/Cp for storage while using local viscosity and conductivity. The original
material planner supplies momentum relaxation and the distinct D3Q7 thermal
relaxation. Native collision families, independent bulk/higher-order choices
and TRT wall-placement parameters are preserved through
`ViscousLattice3::with_shear_relaxation`.

Opening coefficients are held over the original step; endpoint queries admit
and retain the next state without a second transport step or nonlinear iteration.
The existing affine law remains the sole buoyancy force owner. This adapter
adds molecular transport, not variable-mass storage, a source-density buoyancy
law, mechanical work heating, turbulence closure or a new temporal-order claim.
Original heat/source/impulse receipts retain their lattice units.

Endpoint material refusal or cancellation precedes the original joint
flow/force/heat/finite-wall publication. A refused step leaves accepted physical
state and the caller's endpoint-background buffer unchanged; workspace scratch
may contain a refused prefix. Source callbacks must be repeatable; their external
effects remain caller-owned. Existing scalar, tensor and SGS entry points retain
the same numerical transaction and publication point.

Two focused regression groups compare explicit local SI-to-lattice rates and
original source samples on D3Q19/D3Q27, then check conducting-wall heat accounts,
late source refusal, final cancellation and exact retry. Their execution is not
claimed by this implementation checkpoint.

### Independently authored bulk viscosity with material SGS

`FlowField::advance_material_momentum_with_bulk_viscosity_and_sgs` supplies
the same predictor with an explicit original molecular bulk field. Molecular
stress is `2 mu dev(D) + zeta div(u) I` and its physical heat is
`2 mu dev(D):dev(D) + zeta div(u)^2`. The added SGS stress remains purely
deviatoric and uses actual cell M/V; its modeled transfer stays separate from
enthalpy. Original molecular shear viscosity alone sets the SGS applicability
limit, so a large bulk coefficient cannot hide a shear-ratio violation.

The original source, reconstruction, nonlinear momentum solve, force/pressure
update and accepted mass-route transport each run once. All-zero closure
coefficients recover the existing explicit-bulk path exactly. The closed
wall/symmetry scope, frozen molecular-state split and independent energy/work
gates remain. Active pressure mouths still require their own boundary-consistent
closure and work partition before they can be admitted.

Three focused regression groups cover independent graded longitudinal decay,
a solenoidal circulation mode whose response is independent of bulk viscosity,
and zero-model/refusal behavior. Source parsing, scoped formatting and scalar
analytic energy checks passed; native execution remains subject to the
dependency build limitation recorded above.

## Adaptive enthalpy and species exchange through reservoirs

`diffuse_species_with_enthalpy_and_boundaries_coupled_adaptive` reuses the original
state-refreshed exponential-midpoint integrator and positive native
reservoir/interior maps. Each half-interval predictor supplies the midpoint
caloric and transport fields; the full corrector starts from the same original
substep inventory. Only accepted correctors contribute physical boundary heat
or species receipts. Preparations, face queries and pair updates from predictors
and rejected resolutions still consume the whole-call resource allowances.

`SpeciesEnthalpyDiffusionModel` has an explicit boundary partial-enthalpy
callback with the original canonical descriptor and input index. Closed models
retain a default refusal for an active unbound boundary law.
`GasSpeciesEnthalpyBoundarySource` implements that contract using its exact
reservoir/physical-face-temperature binding. Each preparation recovers original
gas temperature from composition and total enthalpy, refreshes admitted
temperature-dependent binary diffusion and transport, and reports its actual
caloric residual. The accepted endpoint source image is returned for reuse.

The independent controls in
`AdaptiveCoupledBoundarySpeciesEnthalpyDiffusionControls` require field accuracy in
K, absolute J/kg and every constituent fraction, plus integrated per-face
accuracy in kg and J. The specific-enthalpy estimate compares both represented
endpoints and retained interval heat divided by carrier mass, so a large
formation datum cannot erase unresolved heat below the stored-h precision. Each boundary compares signed and both gross transfer
histories; the maximum of those three coarse/fine differences is summed over
faces and divided by three. Thus opposing reservoirs and heat reversals cannot
hide unresolved exchange behind a stationary endpoint or zero net domain heat.
These remain smooth second-order estimates, not certified model-error bounds.

Original low components survive both energy and constituent history
accumulation across the complete interval. The initial-to-final mass and
enthalpy audits subtract those retained external transfers before rounding;
formation-reference scale does not relax the absolute joule allowance.
Returned signed and positive/negative NET-map heat is already in the endpoint
enthalpy. Face conductances report the last accepted corrector midpoint,
whereas transfer histories cover the entire interval.

Two independent stiffness checks prevent false coarse/fine agreement.
Changing sampled diffusivity or relative partial enthalpies must resolve the
full OUTER face relaxation to at most one half, regardless of inner sweep
count. When reservoirs share a cell or compete with active interior maps,
the full INNER relaxation must also be at most one half, even for constant
coefficients: otherwise the first saturated map can incorrectly claim all
exchange. The topology check is linear in authored boundary count and allocates
nothing. Isolated exact maps and explicitly resolved fixed-law inner sweeps
retain their stiff paths. The existing frozen-coefficient adaptive reservoir
entry uses the same allocation-free competition check and its native inner
refinement controller. This guard resolves sampled saturation, not an
unsampled constitutive discontinuity.

Boundary histories, native scratch, final readbacks and at most three live
source preparations are covered by checked storage allowances before their
allocations. Cancellation and any source, mass, caloric, work, storage or
refinement refusal publish no partial physical result. The caller continues
to own external reservoir storage, the supplied physical interface temperature
and any pressure/volume evolution.

## Linear characteristic outlets on native D3Q19/D3Q27 flow

`OpenBoundaryGrid3::set_characteristic_outlet` attaches one held incoming
isothermal acoustic characteristic to an existing pressure face. Its uniform
density is the fixed reference; `c = sqrt(L::CS2)` remains lattice-native.
For the outward normal, `J± = (u_n-u_ref) ± c*(rho-rho_ref)/rho_ref`.
The original streamed donor supplies J+ and both tangential velocities;
the control supplies J-. Density and normal velocity can therefore carry an
outgoing pulse together instead of fixing pressure independently of velocity.
Reference, donor and resulting normal speeds must be strictly outward and
subsonic. The caller also bounds relative density departure from the reference.

The same native regularized stress reconstruction, force-density/actual-density
acceleration endpoints, population source frame and measured reservoir exchange
remain authoritative. Serial and pooled paths share the trace map. Material
rheology and SGS use the actual one-sided velocity gradient at these outlets;
ordinary pressure faces retain their zero-normal-gradient condition. Existing
coupled thermal steps consume the resulting candidate physical velocity and
retain their final joint publication point. The outlet adds no step allocation,
callback or retained physical history.

Pressure profiles and characteristic controls are mutually exclusive. Clear
the control before changing its uniform reference density; unrelated inlet or
wall control updates retain it. Mixed wall/open nodes are refused because their
separate population ownership cannot impose both acoustic moments. Backflow,
sonic flow, exceeded density/Mach bounds and cancellation preserve accepted
populations, source epochs, controls and time. Clearing the control restores the
original pressure operator. This is linear, local plane-normal acoustics, not
nonlinear NSCBC, oblique transparency or a compressible caloric gas boundary.
See Poinsot and Lele, [JCP 101 (1992), 104–129](https://doi.org/10.1016/0021-9991(92)90046-2),
for the characteristic framework; no general reflection bound is claimed.

Focused native cases exercise all face normals and collision families, physical
force/stress moments, actual-density acceleration and atomic refusal, reference
equilibrium, and the shared gradient. A traveling acoustic pulse compares the
returning characteristic against the original fixed-pressure termination on
both stencils and outward-normal signs.

### Native source-resolved molecular and SGS convection

`BoussinesqOpenGrid3::{step_serial_with_material_transport_and_sgs,advance_serial_with_material_transport_and_sgs}` compose the original material and local SGS owners in one native transaction. Each accepted opening temperature resolves the original source once; its bounded density/Cp projection and retained SI scaling produce fresh molecular viscosity/conductivity rates. The existing Smagorinsky/WALE/Vreman preparation interprets the opening populations with their own accepted spatial total force, adds local eddy viscosity and `alpha_t=nu_t/Pr_t`, and retains the original collision family, higher-mode rates and TRT wall-placement rule. The material workspace remains molecular; repeated steps never recycle augmented SGS fields as molecular input. The separate SGS workspace carries the actual effective fields.

The original native thermal/force stage advances once. Material queries at actual final temperatures, all source/storage/closure/thermal/force admissions and the common final cancellation gate precede joint publication of flow, heat, finite-wall history, force epochs and receipts. A changing endpoint background remains caller-owned on refusal and exchanges with the old background only on success. Returned evidence contains both concrete material-source epochs, the original heat/force receipt, and the separate modeled-stress and turbulent-heat diagnostics. Both reusable workspaces may contain refused scratch; neither is accepted physical state. Existing material sample/output caps and SGS workspace/ratio caps remain in force. Molecular and modeled coefficients are held over the interval; reference heat storage, affine buoyancy, thermal splitting and the prohibition on depositing modeled or body-force work as heat are unchanged.

## Material SGS through pressure-reservoir mouths

`FlowField::advance_material_momentum_with_viscosity_sgs_and_pressure_profile`
composes original-source molecular momentum, resolved-gradient SGS and the
existing spatial pressure/material-transport owner. Exactly one opposed pair
of reservoir mouths prescribes normal traction and zero tangential velocity.
Actual cell M/V supplies SGS dynamic viscosity and every material-ratio cap;
the original physical boundary half-cell masses still supply inertia.

Cell gradients use the admitted tangential trace. Each mouth shear half-edge
evaluates the selected Smagorinsky, WALE or Vreman closure from its native
reciprocal differences: boundary normal velocity along the mouth and tangent
velocity against the zero trace over the normal half-cell distance. Remaining
tensor entries and filter widths use the original graded-volume interpolation
policy, not a claim of an exact surface gradient. Separate low/high caches
preserve different closures even with one cell between the reservoirs.
Surface-sized storage is checked and fallibly allocated; refresh/admission
polls cancellation. Closure evaluation occurs once per native refresh and
both reciprocal momentum rows reuse the same coefficient.

The two adjacent total viscosities `mu + (M/V) nu_t` enter the original
width-weighted harmonic mouth stress. Its physical half-edge power is divided
by total material resistance, then into local molecular heat and modeled SGS
transfer. Only molecular heat enters h. Every adjacent material independently
constrains the mouth closure, including endpoint molecular admission. The
following original force/projection stage applies the spatial reservoir load
once; pressure work, time/scatter losses and SGS transfer remain separate.

Zero coefficients recover the original molecular reservoir advance exactly.
The entry retains the existing first-order frozen-source split, exploration
policy, atomic refusal and fixed-geometry boundary scope. It supplies no inlet
turbulence spectrum, subgrid-energy storage/cascade equation, wall function or
variable-density total-energy qualification.

Ten focused regression groups are authored: graded boundary traces, all six
normal/tangent permutations and three closures, actual-density edge limits,
heterogeneous traction/heat partition, reciprocal derivative cancellation,
native half-cell decay, bidirectional inventory transport, nonuniform pressure
loading and datum invariance, disabled models and refusal/replay. Pinned
nightly syntax checks, scoped formatting and independent scalar energy
arithmetic passed. Native compilation and execution remain unverified after
the upstream asupersync build limitation; no test-pass or performance claim is
made for these unexecuted regressions.

### Independent molecular bulk response at SGS reservoirs

`FlowField::advance_material_momentum_with_bulk_viscosity_sgs_and_pressure_profile`
adds the authored molecular bulk field to the same original material-reservoir
composition. Molecular stress remains `2 mu dev(D) + zeta div(u) I` and its
physical heat retains the separate `zeta div(u)^2` contribution. SGS remains
purely deviatoric, uses actual M/V, and contributes only to its modeled-transfer
account. Molecular shear viscosity alone constrains the SGS ratio; increasing
bulk viscosity cannot relax that cap.

Normal-traction mouths retain their actual physical half-cell inertia and
shared two-cell shear closure. The source, predictor, pressure loading, heat
insertion and accepted mass-route transport retain their original owners.
Zero coefficients exactly recover the existing molecular-bulk pressure path.

Three regression groups are authored: the independent nonlinear longitudinal
mouth root for three bulk responses, bulk-independent solenoidal mouth shear,
and all zero models plus invalid bulk/cap/refusal replay. Original molecular
heat alone reaches enthalpy; native and original-cell energy identities include
their distinct BE/scatter losses. Pinned syntax/scoped-format and independent
scalar balance checks passed; native compilation/execution remains unverified
under the dependency build limitation described above.

### Source-bound gas Fourier conduction

With `enthalpy-coupling`, `MaterialAdvectionStep::conduct_thermochemical_species`
binds every accepted fluid-cell mass and composition to its original
`GasEnthalpyLaw`, then advances total m*h through the existing native
tensor-caloric backward-Euler or SDIRK2 owner. The original safeguarded
NASA inverse supplies T(h,Y) and the same query's analytic Cp throughout
Newton iteration. Formation energy remains in h; Fourier heat is deposited
once through the native stored-energy update. Zero species diffusivity does
not suppress this independent heat mechanism.

Gas conductivity follows the explicit `ThermochemicalConductivityModel` in
the heat controls. `HeldAtOpening` preserves the opening isotropic source value
and its first-order freezing error. `SourceIntegral` prepares continuously
supported original species viscosity/conductivity paths and integrates their
actual Wilke/WMS mixture throughout the nonlinear solve. The original endpoint
k supplies the native Kirchhoff tangent, without a second dk/dT term.
Endpoint preparation refreshes temperature, conductivity and diffusivity for
continuation. Both SDIRK stages use the same authored thermal setup. Original
heat-flux, face-temperature, film, source, contact and sink semantics apply.
Fixed cell temperatures derive prescribed h from the same gas law and must
satisfy the native inverse-consistency tolerance.

One immutable caloric model is constructed per cell per interval. Trial
queries reuse its original source law without cloning material tables or
resolving transport. Fallibly reserved arrays and a shared caloric-query
allowance bound adapter allocation failure and inverse work; the complete
per-cell model/card storage and native spatial solver storage remain additional.
Every inverse retains the original root iteration ceiling and cancellation gate.
The absolute native interval energy closure plus initial and endpoint caloric
residuals must fit the declared joule allowance. Independently, the sum of each
cell's largest actual trial inverse residual times its mass must fit that
allowance. These source residuals are diagnostics, never extra physical heat.
Local equations retain the native thermal solver tolerances.

This entry admits only fluid voxels with positive accepted carrier mass.
The native heterogeneous caloric API remains the owner for positive-mass
solid storage. Source pressure, composition and carrier inventory are fixed
during this heat component; EOS density is diagnostic. This is not a sealed
rigid-volume internal-energy/Cv pressure evolution. Any source, numerical,
work, energy or cancellation refusal returns no partial candidate.

Original-source regressions are authored for graded unequal-mass BE/SDIRK
heat modes, nonlinear NASA storage and formation-reference covariance,
uniform/source/temperature-boundary limits, fixed-cell storage heat, and
query/source-domain/cancellation refusals. Source parsing and scoped format
checks passed; native compilation and integration execution remain unverified
because the available filesystem is full and DSR/RCH are unavailable.

### Source-bound gas stream advection

`GasSpeciesEnthalpySource::advection_inlet` evaluates canonical mass fractions
and a physical temperature once through the original full source bank and
pressure context. Its private immutable handle binds the resulting density,
bulk total h and composition together. Each boundary/backflow face and local
injection retains its own explicit state; extraction uses the original donor.

`FlowField::project_and_advect_thermochemical_species` admits the original
conserved M, constituent masses and total h, checks initial source closure,
then runs the existing projected material transport once. Constituents consume
that accepted carrier's SAME retained mass routes and operator. Endpoint
temperature and transport are reconstructed through the original source from
the returned inventories; EOS density is diagnostic and cannot replace M.
No stream heat, formation heat or pressure-flow work is deposited again.
The source pressure remains the supplied constitutive coordinate and the
projection's mechanical account remains separate.

The entry retains the native mass, energy, partition, iteration and storage
controls. Additional checked budgets cover lowering arrays, the endpoint
source field and initial/endpoint caloric inversion counts. Original material
evaluation scratch remains additional, and prebuilt inlet evaluation is a
separate source operation. Positive constituents may not disappear through
unrepresentable mass fractions or molar amounts. Initial and final original
caloric residuals plus the native transport defect must fit one absolute
joule allowance, independent of formation-energy scale. Any source, resource,
conservation or cancellation refusal returns no partial candidate and changes
neither the supplied inventories nor the flow.

Four regression groups are authored for independent boundary/local-source
mixing under three formation references, exact native material/species parity,
reverse-face inflow and continuation, immutable source/context binding,
resource/partition/cancellation refusal, exact absence and solid entries.
Pinned source parsing and scoped formatting passed; native type checking and
execution remain unverified under the build/storage limitation above.

## Sensor-protected dynamic WALE

`PeriodicDynamicSgs::with_wale` and `ChannelDynamicSgs::with_wale` explicitly
select the Baya Toda--Truffin--Nicoud shear/vortex-sensor formulation
(ECCOMAS CFD 2010, sections 5.1--5.2, equations 24--26;
https://imag.umontpellier.fr/~nicoud/PDF/ECCOMAS_BayaToda.pdf).
The existing full-gradient WALE kernel replaces the momentum rate in the
original Favre--Germano identity. For `D = dev(sym(G G))` and `S = sym(G)`,
the sensor is `(D:D)^(3/2) / ((D:D)^(3/2) + (S:S)^3)`. Above the authored
threshold, the existing spatial or homogeneous-plane least-squares fit
supplies C_w squared. At or below it, an authored nonnegative fixed C_w
squared is selected before dividing the ill-conditioned near-wall
contractions. Threshold and fixed C_w are explicit; the paper's examples
are 0.09 and 0.5. There is no unprotected dynamic-WALE default or new
coefficient clipping. Existing no-backscatter and coefficient/molecular-ratio
admission policies remain authoritative.

The independent scalar identity retains its original strain rate and
arithmetic. Momentum WALE does not silently change heat or species mixing.
The extra numerical payload is two scalars per cell (WALE rate and SVS):
16 bytes/cell, giving 152 nominal bytes/cell for momentum. Original
Smagorinsky storage remains 136 bytes/cell, with no optional-array allocation.
Actual retained capacities, geometry, finite arithmetic and cancellation
continue through the existing bounded owners.

`NsfGrid3::with_dynamic_wale_thermal_sgs` and
`with_channel_dynamic_wale_thermal_sgs` attach that momentum choice alongside
the original independent heat fit. The required numerical payload is 216
nominal bytes/cell; `required_dynamic_wale_sgs_workspace_bytes` reports it.
`dynamic_wale_config` exposes the selected policy. The periodic/channel
builders preserve the original same-state force correction, current wall
data, half-cell derivatives, conservative energy faces and publication
transaction. No separate PDE, physical timestep or thermal wall law is added.

All 100 actual-source SGS tests passed in a focused native build with the
complete production math crate and original model, Vreman, tangent and dynamic
subtrees. Six new groups use independent physical-space filter/tensor oracles,
check all channel normals and rectangular geometry, the cubic near-wall rate,
unit and Galilean transformations, bitwise preservation of the independent
heat fit, and cancellation/retry with reusable storage. Three NSF owner
regression groups are authored and format-parsed but were not compiled or
executed in this environment; the shared numerical tests do not establish
full NSF/runtime integration. This remains exploratory LES, without a wall
function, shock closure, complete compressible SGS-energy model, DNS
calibration or acceptance-render qualification.

## Nonlinear material rheology with SGS and thermal-volume transport

`FlowField::advance_piecewise_thermal_volume_momentum_with_endpoint_rheology_and_sgs`
extends the existing liquid material/phase owner with SGS at the SAME native
viscous endpoint. Every coupled sweep observes the original unit molecular
quadrature `sqrt(2 dev(D):dev(D))`, resolves the original source's shear and
independent bulk response, then refreshes cell, interior-edge and reservoir-edge
SGS from that velocity before the common residual. The original time level,
actual M/V density and face masses are retained throughout one bounded solve.

Physical rate observation temporarily moves out both the cell eddy field and
the state containing all edge closures, swaps in unit molecular coefficients,
then restores the original fields on success, observation refusal or cancellation.
This retains cache allocations and prevents modeled coefficients from biasing
source operating coordinates, including pressure-mouth half-edges. The first
exact source coefficients are installed before SGS initialization and cap
admission; an arbitrary initialization viscosity never sets applicability.
All later caps use freshly resolved molecular shear viscosity. SGS adds no
bulk response and shares no inferred material law.

The accepted molecular and SGS stresses share one native momentum/work audit.
Only original molecular shear/bulk heat joins authored heat in the existing
caloric owner. Compatible volume change, pressure work and accepted material
transport then run once. Returned physical source rates, original phase
receipts, bulk coefficients, nonlinear work counts and separate modeled-energy
accounts describe that same viscous endpoint. The original temperature is held
during this constitutive solve; simultaneous temperature/rate/SGS feedback is
not supplied by this entry.

The existing fully liquid phase support, source validity, fixed geometry and
optional opposed normal-traction mouths remain. Zero closure coefficients
recover the original endpoint-rheology result exactly after configuration
admission. Source, cap, phase, nonlinear-work or cancellation refusal publishes
no partial material endpoint.

Five focused regression groups are authored: pure-rate observation and cache
restoration; true initial-source ownership of the cap; three independently
derived nonlinear source/SGS roots with original material receipts and
molecular-only calorics; all disabled models; and late source, cap, phase,
work and cancellation refusals with retry. Independent scalar energy identities
were checked. The local executor became unavailable during this increment;
the final source was reviewed and prepared through GitHub, without a completed
Rust formatting, type-check or test-execution gate. The earlier dependency
build limitation also remains. No runtime or turbulent-accuracy claim is made
from these authored but unexecuted regressions.

## Implicit temperature/rate rheology with SGS momentum

`advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs`
composes the existing implicit fluid/finite-solid thermal feedback with native
endpoint-rate molecular shear, independent bulk response and the selected
Smagorinsky, WALE or Vreman closure. Every predictor and thermal trial resolves
the original material source at that trial's exact phase and physical shear
rate; accepted SGS stresses and local molecular-relative caps use those same
coefficients and actual original M/V. Optional opposed reservoir-mouth stresses
retain their native boundary owner.

All trials retain original fluid M/P/h and finite-solid storage. Only physical
molecular dissipation supplies heat; modeled SGS transfer stays separately
reported. Conductivity in this entry remains molecular. The original shared
non-cancelling heat/energy gates and whole-call source, stress and thermal work
ceilings apply across every rejected trial. Force, pressure and material
transport run once after coupled acceptance. Zero coefficients retain the
original thermorheological result; source/phase/cap/work/cancellation refusals
return no partial material or solid endpoint. SGS remains explicit exploration.

Focused regressions use the original temperature/rate material-card fixture and
an independent scalar momentum root with analytically eliminated thermal rows.
They cover joint source coordinates, molecular-only heating, separate modeled
transfer, original finite-wall history, all three zero closures, refreshed
thermal molecular caps, bounded work, cancellation and deterministic retry.

### Original gas and finite solid calorics in one heat solve

`MaterialAdvectionStep::conduct_thermochemical_species_with_solids` assembles
accepted gas inventories and explicit `ThermochemicalSolidConductionCell`
bindings into the existing heterogeneous tensor-caloric BE/SDIRK2 owner.
Every solid voxel has exactly one binding in increasing cell order, with its
own positive reference mass, previous specific energy, original caloric law
and tensor conductivity law. The unchanged gas carrier retains zero solid
mass. The native material table admits occupancy and material-pair contacts;
it does not reconstruct either gas or solid thermal inventory.

Gas and solids share the same physical face equations. Internal interface
heat enters their storage equations once with opposite signs, including
finite face resistance and original source-contact laws. Each solid retains
its declared storage basis, source domain, reference constraints and nonlinear
conductivity integral. Original phase curves retain exact zero dT/de on a
latent plateau; no apparent-Cp smoothing or additional latent deposit is added.
Fixed solid temperatures require explicit prescribed energy because temperature
alone cannot choose a phase fraction. Native fixed-state coverage and original
T(e) consistency checks remain authoritative.

The untouched native result contains physical temperatures and original
specific energies for every voxel, both stages where requested, and the
combined gas-plus-solid interval ledger. A separate gas h vector masks solid
slots to exact zero and feeds original gas endpoint preparation. Initial
combined-energy storage is reused for this masked output. The existing
gas-only entry retains its borrowed endpoint getter.

The shared caloric-query allowance counts gas-source inversions; solid queries
retain original law admissions and native shared Newton/backtracking/Krylov
limits. Gas caloric residuals accompany the combined native energy closure
within the absolute joule allowance. Bindings, source refusal, work limits,
energy failure and cancellation return no partial gas or solid state.
Occupancy, reference masses and gas composition remain fixed during the
interval. Crossing an admitted latent chart does not infer moving walls,
expansion, mass transfer or a new phase geometry. Gas conductivity follows the
same explicit held/source-integral choice as the gas-only entry; each supplied
solid law keeps its original k(T).

Focused original-path regressions cover graded unequal-mass gas/solid heat
exchange, finite interface resistance, latent storage and crossing, prescribed
plateau energy, original solid mass binding and transactional refusals.
Pinned syntax/scoped-format checks passed. Native integration compilation and
execution remain unverified under the previously documented runner/storage
limitation.

### Original nonlinear gas conductivity through the native heat owner

`ThermochemicalConductivityModel::SourceIntegral` supplies the native
`TensorConductivityLaw` from the same cell's already-bound original caloric
model and source context. Prepared pure-species viscosity and conductivity
lines retain every original knot, selected claim, material identity and fixed
auxiliary condition. They are constructed once per cell, outside Newton.
Every point reuses the existing original Wilke interaction and WMS summation;
neither mixture-k interpolation nor a midpoint sample substitutes for an
interval integral. Both gas-only and mixed gas/solid entries use this path.

The material owner splits at every crossed source knot and evaluates bounded
embedded Gauss-7/Kronrod-15 means in normalized coordinates. Its absolute and
relative tolerances have conductivity units; the returned error is a numerical
estimate, not an enclosure or source uncertainty. The estimate does not enter
the physical joule account as heat. Complete original source support is required
over the declared caloric bracket, including interiors; gaps, nonpositive paths
and changing source selections refuse rather than extrapolate.

One point-evaluation allowance spans every fluid cell, face integral, Newton
trial and stage. Preparation uses separate cancellation checkpoints around
bounded original line queries. Each mean retains a shared subdivision limit
across its source segments and a fixed depth/stack cap. The report includes
actual conductivity evaluations and the largest estimated absolute mean error
over accepted and discarded queries. Native nonlinear/linear tolerances,
gas caloric-query limits and physical energy admission remain independent.
The original WMS interaction retains O(N^2) species work and its normal small
receipts per point; no performance benchmark or allocation-free query is claimed.

The source-integral choice removes the held-k approximation within this
fixed-composition, prescribed-pressure model. Smooth SDIRK2 temporal behavior
still requires resolving caloric, spatial, nonlinear and quadrature errors;
unresolved source or phase corners require refinement. Five tests of the actual
production quadrature leaf compiled and passed. Four original-material and
three native spatial regression groups are authored, covering original point
parity, nonlinear means and source gaps, analytic graded heat with crossed
conductivity knots, formation-reference covariance, constant-k BE/SDIR parity
and transactional work/source/cancellation refusal. Pinned parsing passed;
the full material and native spatial tests remain unexecuted.

### Spatial thermal reservoir profiles

`ThermalOpenField3::{set_reservoir_temperature_profile,clear_reservoir_temperature_profile,reservoir_temperature_profile}`
hold independent temperature values on each original reservoir plane. Face
order is XMin/XMax/YMin/YMax/ZMin/ZMax; plane entries are x-fast with the normal
coordinate omitted, matching native momentum profiles. Each incoming fluid
link uses its own authored temperature in the original Danckwerts total-flux
condition. Outflow retains its existing zero normal conductive flux and uses
the actual outgoing fluid state; flow reversal selects that same link's held
reservoir value. Original scalar-population heat, inflow/outflow-volume and
resolved-trace receipts remain authoritative. There is no additional heat
source, thermal state or per-step allocation.

All authored entries, including currently outflow and inactive solid slots,
must pass the original finite scalar-temperature bounds before any input
changes. Native theta may be signed; material adapters retain their physical
Kelvin admission. `BoussinesqOpenGrid3` additionally admits every profile
temperature against its affine density-departure envelope, both when setting
a profile and when attaching a pre-profiled thermal field. Actual candidate
boundary traces still pass the original gates before joint publication.

Profiles share `ThermalOpenLimits3.maximum_storage_bytes` with the original
thermal arrays and other installed profiles. Initial installation allocates
one face vector; subsequent same-face updates reuse its storage. Failed
admission or allocation preserves the previous complete profile and physical
state. The scalar boundary declaration remains an admitted fallback:
uniform-target updates refuse while a profile is installed, and clearing
restores that unchanged fallback without advancing heat, finite walls, flow,
forces, time or history. Thermal face kinds remain fixed; profile operations
on non-reservoir faces refuse. A later refused or cancelled coupled step
retains the previously authored profile and all accepted physical state.

## Native porous through-flow and prescribed skeleton velocity

`lattice3::open_boundaries::forcing::resistance::ResistanceOpenGrid3`
couples the existing implicit `DarcyForchheimer3` law to the original native
D3Q19/D3Q27 velocity/pressure reservoir operator. The held law is
`F = b - D(u-v_s) - beta |u-v_s| (u-v_s)`. Fluid velocity, prescribed skeleton
velocity and material axes share grid coordinates. Coefficients are lattice
force-density/velocity and force-density/velocity², with any physical density
factor already included; the owner does not multiply drag by rho again.

Attachment changes the native half-source representation while preserving
physical density, velocity, force-corrected stress, time and boundary controls.
Collision uses the accepted force. Streamed donors and final populations use
their own density, raw momentum and the original implicit constitutive solve.
A complete reservoir trace computes its law-derived force at the prescribed
physical velocity, attaches the original collision-specific source, and
re-admits the represented endpoint. Spatial/time-updated profiles and the
linear characteristic outlet continue to own their original target moments.

Active linear or nonlinear resistance on a mixed wall/open cell refuses.
Zero-drag seams retain the original force-density reconstruction. Interior
porous cells and stationary voxel walls are supported; exterior motion follows
the existing tangential, locally impermeable fixed-wall law. Rebinding controls
cannot change topology or expose a mutable force-only step that loses drag
history. Constitutive coefficients and skeleton motion remain held.

The original `ResistanceMomentumExchange3` reports independent body impulse,
skeleton reaction/work and nonnegative relative dissipation, plus actual
reservoir transport, ALL stationary/moving/voxel wall links, and native periodic
angular transport. Pure reservoir links are excluded from the wall probe;
mixed links retain wall ownership. No term is inferred from a global residual.
Drag is retained separately from total force so a larger independent body load
does not erase its constitutive receipt. Original unresolved half-force and
work-cancellation refusals remain authoritative.

Native stream tiles resolve force and drag into reusable scratch; serial
steps allocate no source storage and invoke no user law callback. A canonical
measurement pass uses accepted/staged physical velocities and the original
periodic pull, without another collision or nonlinear solve. Populations,
force/drag history, time and the last complete receipt publish only after all
endpoint, wall, work and cancellation checks pass. Existing source-density and
actual-density-acceleration paths retain their original arithmetic.

This is fixed-geometry, low-Mach isothermal porous flow with prescribed
skeleton motion. It does not evolve porosity, remap moving solids, integrate
skeleton inertia or deposit drag heat. The held-law work identity
`fluid_drag_work + skeleton_work = -relative_dissipation` is not a total
transported fluid kinetic-energy balance. The focused tests cover preserved
attachment, a manufactured anisotropic/nonlinear physical endpoint and common
velocity boost, pressure-driven porous-slab flow, zero-law population parity,
complete wall/reservoir/periodic momentum, updated traces, seam refusal and
transaction rollback.

## SGS heat transport at the coupled material endpoint

`advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs_transport`
adds explicit turbulent-Prandtl heat transport to the same original-history
implicit temperature/rate/SGS solve. The caller supplies finite positive Pr_t.
Each fluid side evaluates k_total(h) = k_source(h) + (M/V) Cp(h) nu_t / Pr_t
from the original material/caloric source and the accepted native cell closure.
Total conductivities enter the original two half-cell resistances before their
harmonic combination. The existing native sparse heat graph owns conservative
fluxes and nonlinear storage; its attached link laws include the complete
enthalpy derivatives of both k_source and sensible Cp. No finite-difference
material derivatives or separate thermal solver are introduced.

The heat graph freezes one closure per trial. The outer owner refreshes momentum,
molecular shear/bulk and SGS at that exact thermal endpoint and charges BOTH the
molecular heat mismatch and twice the sum of absolute modeled face-heat changes.
Both adjacent residuals are charged, so opposing face exchanges cannot cancel.
The same bound enters the explicit feedback tolerance and final shared thermal
budget; the independent mechanical budget still accounts for molecular heating
and modeled momentum transfer through its original owner. Rejected trials relax
heat and closure while retaining the original fluid and finite-solid histories.

The feedback report retains Pr_t, the closure used by the accepted heat graph,
the accepted native closure and the explicit modeled heat-coupling defect. Native
coefficients transfer by ownership after the successful stress/work audit; the
thermal graph does not infer them from scalar dissipation. Validated zero model
coefficients take the original path with empty closure receipts.

This isotropic gradient-diffusion closure applies only to fluid-fluid links.
Baths, fluid-solid interfaces and source contacts retain their molecular laws;
no thermal wall function, scalar-variance equation, variable-Pr_t model or
immediate conversion of modeled SGS energy transfer into caloric heat is claimed.
Positive modeled transport requires a supported fully liquid sensible response.
Original material domains, exploratory SGS policy, cancellation and whole-call
thermal/stress work ceilings remain authoritative.

Focused regressions cover the genuine variable material response, total
half-cell resistance and analytic tangents; independent coupled momentum/thermal
roots at two Pr_t values; finite-solid history and molecular wall interfaces;
used/accepted closure receipts; exact zero closures; invalid Pr_t; bounded work;
late cancellation and original-input replay.

### Radiation on original gas/solid caloric storage

With `enthalpy-coupling`, `advance_caloric_tensor_conduction_radiating`
and `advance_caloric_tensor_conduction_radiating_sdirk2` couple the original
backward-Euler and SDIRK2 tensor/caloric owners to the existing gray
radiation geometry and checked enclosure radiosity. Every radiative retry
retains the same physical starting specific energies and reference masses.
Original nonlinear conductivity, caloric slopes and source contacts remain
live; an isothermal latent plateau changes stored energy without an extra
latent heat deposit. Original Newton, Krylov, backtracking and source query
allowances are shared across all retries, and the explicit radiation budget
bounds complete native heat solves.

Both the temperature update and original physical radiative power must
converge. The native endpoint retains the final tangent-sink diagnostics;
`TensorCaloricRadiationReport` discloses their physical-power discrepancy.
The interval record uses the actual emission/absorption correction.
Prescribed-cell radiation changes the thermostat reaction by the same amount,
so it does not introduce an additional net external heat term. Any invalid
source, inventory, work budget, unresolved radiation or cancellation returns
no accepted interval.

`MaterialAdvectionStep::conduct_thermochemical_species_with_solids_radiating`
uses the same mixed original gas/solid preparation, original species composition,
finite masses, and gas caloric/source-conductivity budgets. Its returned
`ThermochemicalConductionAdvance::RadiatingBackwardEuler` or
`ThermochemicalConductionAdvance::RadiatingSdirk2` retains both native heat
and radiation evidence; `radiation_reports()` borrows the stage reports.
Gas-only continuation still uses the separately masked original gas enthalpy
field. The selected backward-Euler or SDIRK2 scheme drives the same original
mixed constitutive objects throughout the interval.

SDIRK2 converges physical radiation separately at both original stage times
using one prepared radiation geometry and one interval work budget. The second
stage carries actual first-stage radiative heat on free cells; prescribed-cell
carry retains its original storage/thermostat convention. The physical heat
record uses the original tableau weights and the same initial energy in both
storage equations. `record.radiated_w` remains the final endpoint power. The
integrated radiative energy is duration times `(1-gamma)*P1 + gamma*P2`, using
the two reports. Smooth-segment second-order accuracy does not apply across an
unresolved phase corner or boundary discontinuity.

The reused radiation discretization emits at exposed solid voxel CELL
temperatures, with transparent gas and the original gray/diffuse/opaque
surface assumptions. It does not resolve emission at tensor face traces,
participating-media absorption, moving geometry or phase-transfer mass.
Radiation and conduction are declared parallel heat routes; callers must not
include the same radiative heat in an authored sink or source. The ray geometry
and nonlinear residuals are numerical estimates, not certified error bounds.

Focused G0/G1/G3/G4 regression definitions exercise independent quartic energy
balance, exact latent plateaus, reciprocal enclosure transfer, prescribed
latent inventory, source/budget/cancellation refusal, and mixed gas/solid heat
with the original source-conductivity option, both SDIRK2 physical-stage roots,
and changing prescribed latent inventory. Changed Rust sources were parsed
with the pinned compiler's formatter. Full native integration execution was
unavailable: DSR and RCH are absent and the local workspace has no build space.

## Original-source Fourier integrals in conservative NSF gas energy

`ThermochemicalNsfGas::with_source_conductivity` explicitly selects
`NsfSourceConductivityControl` for the existing NSF total-energy owner.
The original midpoint conductivity remains the default. Source mode retains
the original selected pure-species viscosity/conductivity lines and WMS
mixture, with an explicit quadrature tolerance and maximum actual WMS
evaluations per face or complete nonlinear boundary root.

All molecular Fourier paths use the selected convention: internal fluid faces,
prescribed reservoir faces, global/local isothermal walls, finite ambient
films, characteristic energy fluxes, and the characteristic incoming-enthalpy
temperature root. A film or characteristic root prepares its source interval
once and shares its actual point-query allowance across every temperature
candidate. An exactly isothermal source face uses one original point query;
it does not construct an artificial epsilon bracket.

For a wide internal face, source mode applies the native temperature stencil
to the original conductivity primitive: `[27*integral(Ti,Tj,k) -
integral(Tm,Tp,k)]/24`, with the native unit conversion and one held face
pressure. It does not multiply a single secant conductivity by the wide
temperature gradient. The original SGS contribution retains its separate
temperature-gradient convention and original molecular/model ratio admission.
Fallback and boundary faces retain their original full-/half-cell lengths,
orientation and single-valued flux ownership.

The held pressure is the existing actual adjacent-pressure interpolation for
that face; walls retain the adjacent gas pressure. Source preparation requires
continuous admitted support over the complete involved temperature bracket.
The native sensible internal-energy chart, density-derived pressure, pressure
force, enthalpy-flow work, stress work and conservative total-energy flux remain
with the same NSF owner. There is no extra heat deposit or change of species
composition. Default midpoint expressions and arithmetic order are preserved;
the selected source policy enters the state identity.

Quadrature error is an embedded conductivity estimate, not a certified energy
bound or additional physical joules. The original face-pressure interpolation,
spatial discretization, fixed-composition and continuum/Mach limitations remain;
this does not establish new general fourth-order accuracy for variable pressure
or experimental flow validation. Existing native cancellation checkpoints
remain, with additional face/root work bounded by the declared evaluation,
quadrature and source-knot limits.

Five G0/G1/G3/G5 native regression groups cover exact piecewise-integral
interior fluxes, finite-film/isothermal equivalence, characteristic upstream
enthalpy counted once, source/work refusal with rollback, and formation-datum
and worker-count replay. All eight changed Rust files were syntax-parsed with
the pinned formatter. Independent scalar integral and film-root calculations
were checked; the native regression groups were not executed because DSR/RCH
and a usable local build graph/space were unavailable.

## Bounded midpoint transport of SGS history

`LagrangianAdvection` explicitly selects original Euler tracing or midpoint
RK2 with a positive authored u8 substep count. The same policy transports
signed stress, independent heat and common-component contraction histories,
including the advected stress used by adaptive momentum memory. Each midpoint
stage samples the current physical velocity through the original positive
trilinear stencil. Channel half-cells interpolate to the actual moving-wall
velocity. Accepted history retains its original positive interpolation and
is remapped only once at the final departure, so trajectory substeps do not
introduce repeated history-remap diffusion.

Every stage and accumulated per-axis travel must satisfy the existing
departure limit, and every sampled stage must remain inside physical walls.
Each characteristic trace is bounded to 255 substeps and at most 510 velocity
evaluations, with cancellation checks at most eight substeps apart. Joint
histories may each trace their departure. Selection precedes the first accepted
sample; trial forks preserve the policy, and rejected work cannot publish
history, clock or normalization. No full-grid arrays are added and original
numerical byte budgets remain unchanged.

NSF exposes explicit policy selection through
`with_lagrangian_thermal_sgs_advection` and reports it through
`lagrangian_sgs_advection`. The original builder retains Euler. Existing
force-corrected current velocity, moving-wall refresh, independent relaxation
times, molecular wall heat, endpoint transactions and reversible native
intervals remain authoritative.

All 106 actual-source shared SGS tests passed in the focused native build
using complete production math and the original model, Vreman, tangent and
dynamic subtrees without numerical or runtime stubs. Six new groups include
independent curved-characteristic and actual emitted-SGS-viscosity refinement,
policy agreement across stress/heat/species, moving-wall interpolation on all
axes, physical-unit changes, unchanged payloads and every-poll cancellation
with exact retry, including adaptive momentum memory. One additional NSF
runtime regression is authored and format-parsed but was not compiled or
executed in this environment.

Second-order convergence applies to smooth frozen-field trajectories only.
Unsteady velocity/source sampling, history-remap accuracy and the complete
physical solver are not upgraded or newly qualified. This adds no inflow
memory, large-Courant qualification, SGS-energy model or DNS calibration.

## Porous momentum and native thermal velocity epochs

`ResistanceOpenGrid3::{prepare_thermal_field, step_serial_with_thermal,
step_serial_with_thermal_tensor_diffusion}` binds the original D3Q7 thermal
carrier to the accepted implicit porous momentum owner. No independent
force-only step or replacement advection/collision kernel is introduced.
The returned `ResistanceOpenThermalStep3` retains the complete original
resistance momentum/work receipt and independent thermal balance.

Opening thermal collision, flux-history correction, reservoir transport and
the first source half-kick use each cell's accepted physical velocity,
including its retained implicit half force. The porous owner stages its
original collision/pull, constitutive endpoint solve, complete wall/periodic
account and work admission once. The closing heat half-kick then uses the
staged populations and their actual resolved endpoint force. Raw population
momentum and a frozen opening force are not substituted for that velocity.

Scalar relaxation and uniform/per-cell SPD diffusivity tensors retain the
existing thermal definitions. Thermal tensors are resolved in grid axes
independently of porous principal axes. Existing reservoir/backflow data,
stationary finite-wall conduction, immersed storage, temperature/kinetic
budgets and moving-finite-wall refusals remain authoritative. Geometry,
solid masks, boundary ownership and accepted step counts must agree.

The composition reuses all native and thermal scratch. Populations,
force/drag history, finite thermal state, flux history, both receipts and
time counters publish only after every momentum, thermal and cancellation
gate succeeds. A thermal refusal after staging momentum preserves both
accepted states and permits the same deterministic retry.

This is passive fixed-reference-capacity heat transport under held porous
laws and fixed geometry. The caller's heat source remains explicit; relative
drag dissipation, skeleton work and mechanical wall work are reported without
an inferred fluid/solid heat partition. Temperature-dependent resistance,
porosity evolution, compressible enthalpy and a total transported fluid-energy
balance are outside this composition. The source ordering is the existing
partitioned native ordering, without a new temporal-order claim.


## Original nonlinear material density in native open-flow buoyancy

`BoussinesqOpenGrid3<L, Q, P>` retains an immutable typed force policy with
its accepted populations, temperature and force epoch. The default remains
`BoussinesqLaw3`, including its existing affine API and admission order.
`BoussinesqDensityLaw3<S>` samples the original `DensityDepartureLaw3`
at actual native temperatures and applies rho_ref times the admitted relative
density departure times gravity, plus the separately authored mechanical
background. Nonlinear and nonmonotone density are supported without fitting
an expansion coefficient. Original source errors remain typed.

The source policy admits each actual source half-kick, transported cell sum,
individual reservoir/external-wall/immersed-wall fluid trace, force preview
and final temperature. Extrema alone cannot certify a nonmonotone relation.
Opening source forces must still agree with retained accepted forces; source
changes cannot silently reframe momentum. The original collision, Guo force,
D3Q7 transport, finite-wall heat, tensor diffusion, impulse account and final
publication gates remain authoritative. No native per-step allocation or
replacement numerical solver is introduced. Refusal preserves accepted
populations, thermal/force/wall history, time and receipts.

`material_density_law3` retains the original reference material owner,
frozen transport snapshot and immutable Kelvin resolver. It uses the original
MaterialScaling SI acceleration conversion and the original Boussinesq
rho/Cp admission, with actual resolved source density in buoyancy. Signed
nonzero Kelvin scales are supported; each actual Kelvin sample must be
positive and finite. Material viscosity/conductivity and material/ordinary
SGS facades support the same typed policy and retain their original source
and reusable-workspace rules. Compositions must supply the same reference,
coordinates and physical source; independent callbacks are not assumed equal.

This remains fixed-reference-mass, fixed-reference-capacity Boussinesq
transport within explicit small-departure limits. Population pressure density
is independent of the material density projection. Variable-mass advection,
compressible total energy, pressure work and a new temporal-order or empirical
fluid qualification are not supplied by this extension.

Eight new regression groups cover density maxima, SI/Kelvin conversion,
nonmonotone boundary and source-step refusal, original material/SGS agreement,
conservation, retained ownership, cancellation and exact retry. All changed
Rust sources passed the pinned formatter/parser. The exact standalone force
law module and its five inline tests passed metadata-only Rust compilation;
the tests were not executed. Full native integration compilation and runtime
validation remain pending in a usable workspace build environment.

## Adaptive original-enthalpy radiation

`advance_caloric_tensor_conduction_radiating_adaptive` composes the original
caloric SDIRK2 radiation stages with the existing embedded energy controller.
It shares the nonradiating adaptive loop, transition schedule, stage forcing,
prescribed-mask admission and whole-window heat gates. The error estimate uses
original specific-energy increments; a sharp latent plateau remains visible
even when its temperature does not change. It is a local asymptotic estimate,
not a certified endpoint or individual heat-channel quadrature bound.

The original ray geometry is prepared once for the entire window and borrowed
by every accepted and rejected attempt. Geometry, emissivities and radiative
surroundings are fixed over that window. Authored thermal forcing may vary at
the actual stage times, with explicit left limits at declared transitions.
Callers needing changing radiative geometry or surroundings must separate
windows at those changes.

`TensorCaloricRadiatingAdaptiveSolution` retains one pair of physical
radiation reports per accepted interval, aligned with `thermal.records`.
Its compensated `radiated_energy_j` uses the native SDIRK stage weights and
is already part of the thermal boundary account; it must not be debited again.
The final record's `radiated_w` remains endpoint power. Completed rejected
trials contribute to work counters but never to accepted heat or history.

Native nonlinear, linear and radiation limits are shared within each attempt;
the attempt limit bounds total work, including incomplete solver failures.
Caller-owned caloric and conductivity query budgets persist across attempts.
Fatal source, invalid-data, allocation and cancellation errors return no
accepted prefix. Whole-window closure and the sum of absolute accepted native
closures must both pass the declared joule tolerance.

G1/G3/G4 regression groups cover nonlinear cooling refinement against an
independent scalar reference, latent energy adaptation at constant temperature,
forcing transitions and prescribed phase inventory, and late work/source/
cancellation refusals. These new native tests are authored and syntax-parsed;
they have not been compiled or executed in this environment.

### Adaptive scheduled evolution of original laminar films

`film::evolution::adaptive` advances the original conservative film owner
through authored continuous drive segments. Each comparison uses one original
backward-Euler interval and its two actual representable halves. Only the
native fine trajectory is accepted; no extrapolated height, replacement
mobility or rejected-trial heat is published. The existing
`ConductionTimeOrder::First` supplies unequal-span error correction and
accepted/rejected step factors through the existing shared thermal dependency.

Local error covers height changes and every original face's signed integrated
volume, viscous heat, body work, surface work and individual reservoir work.
Channels are compared before summation, so constant height and opposing work
do not hide their sampled transfer discrepancy. These are local indicators,
not rigorous global bounds or proof that arbitrary source oscillations,
spatial error or constitutive uncertainty have been resolved.

The drive query carries an exact time and segment index. A query at a segment's
closing event selects its left limit; the next interval selects the following
segment. No comparison crosses an authored event. An endpoint exactly one
floating-point successor below an event or target snaps to that exact time,
using its actual span in all native solves, work and error calculations.
`suggested_dt_s` exposes the final comparison-width proposal for continuation.

Every Carreau query explicitly supplies temperatures for all original profile
sources, including both open ends. Rebinding retains each original law,
support and shear envelope; identical temperature bits reuse preparation.
Newtonian source owners require no temperature vector and retain their frozen
properties. Density, surface tension and normal gravity remain fixed. This is
prescribed source-temperature evolution, not thermal storage, an EOS or
thermal expansion. The quasi-steady, fully wet, small-slope film and its
original laminar-regime admission remain unchanged; fast momentum transients,
contact lines and rupture require other owners.

Native endpoint/profile ceilings apply to the complete adaptive call,
including coarse, rejected and failed work. Only temporal inaccuracy and
original Newton nonconvergence may refine. Source, physical regime,
conservation, cancellation and work refusals remain terminal. Final volume and
stored-energy changes are independently recomputed from opening/final heights
and reconciled with accepted-only transport/work. Original backward-Euler
numerical defects remain separate from viscous heat. Any refusal preserves
heights, source temperatures and time, including after a private accepted
prefix. The material wrapper retains the same borrowed original source and
receipts while using this native adaptive path.

Four native regression groups are authored for analytic throughflow and
counterflow transport/heat, gravity/capillary mode decay, exact decimal events
with original Carreau temperatures, and work/source/cancellation rollback.
A source-bound integration regression covers both original viscosity bases,
native parity and retained source identity. After executor recovery, the
existing original-source recipe compiled the current film, material, math,
solver and shared time-order sources without warnings. All 13 evolution tests
passed: four new adaptive groups and nine existing groups, with 33 unrelated
tests filtered. Rust formatting and focused whitespace checks passed. The
separate full material-source wrapper integration remains authored and
unexecuted. This is focused original-source verification, not a full Cargo,
DSR or asynchronous-runtime integration result.


## Spatial force epochs and reservoir targets in one accepted step

`ForcedOpenBoundaryGrid3::{step_serial_with_open_boundaries,
advance_serial_with_open_boundaries}` prepares candidate reservoir targets
through the original open-boundary owner and stages the original force-aware
collision/stream/reconstruction. The held variant borrows the accepted force
field; the changing-source variant uses the accepted opening force and supplied
endpoint field at their original cells. Targets publish only after native
transport/wall measurement and final cancellation admission, alongside
populations, force history and step count. Refusal retains all accepted epochs.
Reservoir transport and trapezoidal body impulse remain separate; no additional
source reframe or physics kernel is introduced. Cross-crate regression
specifications are in fs-couple's scenario_open_boundary target. This addition
passed pinned Rust parsing and diff checks only; typechecking and runtime
execution remain unverified.

## Adaptive source-bound gas and solid enthalpy exchange

`MaterialAdvectionStep::conduct_thermochemical_species_with_solids_adaptive`
advances the existing mixed gas/solid heat owner to an exact supplied target
with `ThermochemicalAdaptiveConductionControls`. Its optional radiation
selects the native adaptive radiating path. Both choices retain the original
caloric laws, explicit carrier and solid masses, conductivity laws, contact
resistances and single-valued interface heat.

Preparation occurs once for the complete window. Original gas caloric and
conductivity query budgets span every accepted and rejected attempt.
`SourceIntegral` queries the original prepared WMS conductivity throughout
that window. `HeldAtOpening` deliberately freezes gas conductivity at the
window opening; temporal refinement does not correct that freezing error.
Composition, source pressure/context, geometry, setup and prescribed energies
remain fixed. The native forcing API owns time-dependent thermal schedules.

The result retains the full mixed thermal endpoint and a separately masked
gas enthalpy field prepared through the original source. Latent solid energy
never enters gas-only slots, and EOS density never replaces carrier mass.
Native `Adaptive` and `RadiatingAdaptive` variants expose accepted records;
`radiation_reports()` borrows their accepted stage reports without copying.
The adapter's native energy defect uses the larger of absolute whole-window
closure and the sum of absolute accepted interval closures, then admits it
together with opening/final source inversion defects. Trial caloric residuals
retain their separate original admission and are not physical heat deposits.

A late source/work/energy/cancellation refusal returns no partial gas or solid
endpoint. Native line search distinguishes an out-of-support candidate from a
fatal original source/work refusal; exhausted source budgets cannot turn into
backtracking or timestep retries. Native solve errors preserve their original
type; other adaptive refusals retain their full diagnostic. The coupled regression compares original
constant-Cp gas cooling against an independent exponential solution while a
solid remains on its latent plateau, with and without radiation, and exercises
shared caloric and conductivity query-budget exhaustion. A native regression
also checks immediate source-budget refusal inside a Newton trial at the
timestep floor. These tests are syntax-parsed but not natively executed in
this environment.

## Spatial and time-updated native wall velocity

`OpenBoundaryGrid3::set_wall_velocity_profile` resolves prescribed tangential
velocity at each actual halfway-link point `p + 1/2 - e_q/2`, including diagonal
links and wall/open intersections. Both D3Q19 and D3Q27 retain the original
moving reflection; every reflected population, relative wall impulse, torque
and work uses that same link sample. Rheology and SGS use the corresponding
straight-normal link at their halfway-wall gradient sample. Spatial-force and
implicit porous owners expose the same controls without detaching force or
drag history. Their scalar/tensor thermal compositions retain physical opening
and closing velocities and the original moving finite-wall refusal.

Sampling occurs only during explicit control updates, in canonical plane-cell
and native-direction order. Two reused face buffers hold Q velocity vectors
per plane cell; inactive solid/non-crossing slots are zero. The caller's byte
cap covers both buffers of all attached wall profiles. Stepping allocates no
profile storage and invokes no source callback. Caller owns units, temporal
sampling and spatial resolution; geometry, interior voxels and skeleton motion
are not changed by a prescribed exterior drive.

Tangency, native collision/grid Mach limits, shared-corner compatibility and
zero moving mass source at each fluid cell must all hold. Tangential fields
with incompatible discrete surface divergence refuse without correction or
reservoir compensation. Normal/swept wall motion still requires a geometry
owner. A failed source query, storage check, physical admission or cancellation
retains the previous complete trace and accepted fluid state. Clearing a trace
re-admits its uniform fallback with the remaining walls. Uniform rebinding on
an actively profiled face refuses; unrelated controls preserve all profiles.
An accepted control update remains held if a subsequent fluid step refuses.

Focused regression definitions cover exact uniform-trace/native parity and
analytic Couette flow, independent spatial first-step velocity and work,
time-updated original closed-wall parity, wall/open link ownership with pooled
replay, physical rheology gradients, active implicit-drag/thermal endpoint
velocities, finite-wall refusal and atomic source/control rollback. Mechanical
wall work remains separate from prescribed sensible heat; this extension adds
no total-energy, high-Mach or moving-geometry accuracy claim.


## Spatial stateless walls on native thermal and buoyant flow

`ThermalOpenField3::set_wall_boundary_field` and the native Boussinesq owner
accept the existing `ThermalWallField` on each external wall plane. Individual
links may prescribe temperature, signed inward heat flux, ambient Robin film
or insulation. The original reflected pair and scalar transport compute actual
surface temperatures and heat; the existing wall vectors are reused without
allocation or a new time-step kernel. All values and shape/ownership constraints
are admitted before mutation. Explicit nonzero spatial flux on an inactive
solid-facing slot refuses. Finite solid storage and conduction ownership cannot
be created or replaced through this input path.

`wall_states` exposes the actual installed field; `boundaries` retains the last
uniform declaration for face ownership. A subsequent uniform setter replaces
the whole field. Binding a buoyancy owner admits actual per-link prescribed
temperatures, and nonlinear source-density policies continue admitting every
actual reflected trace. Accepted populations, thermal/force/wall history,
receipts and clocks retain their original transactional publication rules.

Four G1/G3/G4 regression groups cover independent reflection/storage on all six
faces, uniform parity, nonlinear interior-temperature refusal, finite/inactive
ownership, cancellation and exact retry. Pinned Rust 2024 formatting/parsing
passed; native integration compilation and runtime execution remain unverified.
No new temporal order, variable mass, moving-solid storage or calibrated
heat-transfer accuracy is claimed.

## Dynamic turbulence at prescribed reservoir faces

`NsfGrid3::with_reservoir_dynamic_thermal_sgs` attaches exploratory
instantaneous Smagorinsky momentum and independently fitted heat transport to
an all-fluid NSF grid with exactly one nonperiodic axis, two complete
prescribed reservoir faces, periodic transverse axes, and at least three
cells per axis. The shared owner retains the original nine-cell transverse
Favre/Germano filter, explicit squared volume-width ratio `9^(2/3)`, and
200 nominal numerical workspace bytes per cell. Domain averaging fits each
transverse plane; LocalBox averaging stays within the transverse neighborhood.
No test-filter or averaging sample crosses an opening.

`sgs3::dynamic::reservoir::{ReservoirDynamicSgs, ReservoirDynamicScalar}`
exposes the same bounded numerical fit independently of the NSF owner.
Velocity and scalar normal derivatives use the current primitive ghosts
one full cell spacing beyond the boundary centers, including normal
throughflow. Boundary replacements refresh the closures through the original
force-rebase and conservative energy-flux transaction. Topology removal
refuses. Cell, shared-face, and prescribed-reservoir heat-ratio admissions use
their original molecular material state; the reservoir check occurs before
extrapolated adjacent turbulent conductivity is added. Original mass,
pressure, energy, trial, and rollback owners remain authoritative.

This scope excludes localized apertures, characteristic boundaries, imposed
unresolved inflow spectra, SGS-energy storage, and qualification as a general
three-dimensional open-boundary LES model. It retains the original passive
eddy-transport restrictions and explicit coefficient/material limits.

Native verification passed five new reservoir groups and twenty existing
channel/scalar groups using unchanged production numerical code and real
`fs_math`. Unrelated test registrations were disabled only in a temporary
snapshot. The initial complete shared-test-tree compilation was killed with
signal 9 under the environment memory limit; the smaller 25-test build
compiled and ran successfully. Two native NSF integration groups specify
current-ghost coupling, force-chart preservation, native-step replay,
cancellation, actual reservoir-face heat-ratio refusal and boundary-update
rollback. Both groups now pass within the 44-group native NSF validation
described below. The full NSF runtime and Cargo workspace were not rebuilt.

Equation and filter-convention reference: Cabot (1991), CTR Annual Research
Briefs, section 2.1,
https://web.stanford.edu/group/ctr/ResBriefs/1991/18_CABOT.pdf .
The homogeneous-direction filtering and independent momentum/scalar
contractions follow this framework; the nine-cell box and volume-width
ratio are explicit discretization choices, not a reproduction of its
spectral filter or empirical validation results.

## Spatial prescribed reservoir faces

`NsfGrid3::set_open_face_profile` installs a complete, already-open reservoir
plane from one admitted `NsfReservoirTileState` per adjacent fluid cell. The
explicit tile ceiling bounds allocation; duplicates, omissions, wrong-face
cells and inadmissible material/Mach/temperature states refuse. Accepted tiles
have canonical transverse order (lowest remaining Cartesian axis fastest),
constant-time lookup and immutable shared snapshots. `open_face_tile_states`
exposes the accepted primitive targets. Installation rebases the original force
chart and energy-flux history while preserving physical density, velocity and
total energy; failure or cancellation retains the previous profile and complete
accepted state. Uniform profiles use the original uniform path.

Normal pressure/energy faces read the local full-spacing ghost primitive.
Incoming diagonal lattice populations read the actual transversely shifted
ghost donor, with periodic transverse wrapping; the receiving-cell correction
is reused and only the requested ghost population is evaluated. Profiled
instantaneous dynamic momentum/heat SGS uses the same original transverse filter
for physical and test-level fields: exterior test primitives receive the matching
nine-cell Favre filter. Existing closure fitting, material and turbulent-
conductivity limits, collision, native mass transport and conservative energy
owners remain authoritative. No extra full-grid SGS halo is allocated. The
topology remains one complete nonperiodic axis with periodic transverse
directions. Prescribed target mass flux is distinct from the compressible mass
actually streamed.

Validation: 28 actual-source shared SGS tests and 44 disjoint actual-source native
NSF tests passed, including independent full-tensor/Favre boundary oracles,
diagonal periodic-donor populations, uniform equivalence, cancellation/retry
and conservative native evolution. The native source build uses real material
and numerical dependencies and omits unused runtime-dependent observation hooks;
it is not a full Cargo/workspace, DSR or cross-ISA gate.

## Authored mass-free dilatation in material and gas transport

`MaterialVolumeSource::Dilatation` and
`GasSpeciesVolumeSource::Dilatation` carry a finite signed integrated volume
rate in m³/s per fluid cell. Positive values require net outward volume flow;
negative values require inward flow. The original pressure projection admits
that constraint, including its closed-domain compatibility and mechanical
work account.

The generic material, force/momentum, resistive momentum and thermochemical
species paths lower the same rate into that existing projection. Dilatation
creates no local transported-mass route, source velocity, source species,
enthalpy injection or additional route-budget charge. Actual corrected
internal and boundary face fluxes retain the original implicit donor mass,
species, momentum and mass-times-enthalpy transport. Boundary inflow still
requires its original incoming physical state; boundary outflow removes
actual donor inventory. Solid cells continue to require `None`.

This is an explicitly authored volume constraint, not an inferred density
change or a new EOS closure. The existing phase-volume owner remains the
authority for its coupled constitutive expansion and volume-admission model.
No projection work or expansion rate is deposited as additional heat.
The original first-order donor and operator-splitting limitations remain.

Focused native regressions cover conservative redistribution under sealed
equal-and-opposite dilatation, actual open-boundary expansion/contraction,
species transport, momentum transport and incompatible sealed net expansion.
These tests are authored and syntax-parsed; native compilation and execution
are unavailable in this environment.

## Spatial current feedback on original caloric storage

`advance_caloric_tensor_electrothermal` connects the original per-cell
electrical/contact solve to the native tensor-caloric backward-Euler owner.
The resistivity callback receives occupancy, original specific energy and
temperature, allowing electrical response to follow phase fraction while
temperature remains fixed on a latent plateau. No new current or thermal
discretization, apparent heat capacity or fitted phase curve is introduced.

Every coupling proposal starts from the same physical opening masses and
energies. Only its proposed Joule sources change. The original conductivity
tensors, physical face traces, contacts, reservoirs and caloric inverse remain
with the native heat owner. A discarded coupling pass adds no physical time or
heat. Endpoint current is reused as the next source estimate, requiring one
more electrical solve than the number of thermal passes.

`CaloricElectrothermalConfig` bounds coupling passes and retains native
per-solve work limits. Acceptance compares energy changes scaled by increments
from the opening state, local Joule-power changes, and an absolute joule
allowance. The allowance includes noncancelling local source discrepancies,
authored/Joule source-combination and reduction rounding, applied/endpoint
electrical work versus deposited heat, and native thermal closure. Existing
double-double arithmetic supports these accounting reductions.

`CaloricElectrothermalStep` distinguishes the electrical field actually
applied from the field refreshed at the accepted thermal endpoint. The latter
deposits no additional heat. Applied Joule heat on prescribed-energy cells is
reported as `thermostat_joule_energy_j`: an equal source and immediate
thermostat-removal pair omitted by the original native thermal ledger. It is
already included in total applied Joule energy. Free-cell Joule heat is already
in the native source account; neither channel is deposited again.

This is stationary, fixed-mass DC coupling with fixed electrode/contact data
and first-order physical time stepping. The bounded source iteration may
refuse strong feedback that does not converge. It adds no expansion work,
advection, chemistry or electrical energy storage. A source, work, energy or
cancellation refusal returns no partial interval.

G1/G3/G4 regressions cover current-dependent latent energy at constant
temperature under an energy-datum shift, spatial single-phase conduction
against a three-cell energy system, contact heat partition with a prescribed
electrode, and exhausted/late-refused coupling. The changed Rust sources are
syntax-parsed and independent scalar references checked; native compilation
and test execution are unavailable in this environment.

### Scheduled adaptive thermomechanical evolution of resolved laminar films

`ThermalCarreauFilm::advance_adaptive` advances the same fixed-inventory native caloric and wall-normal Carreau momentum owners. Each comparison runs one coarse interval and two actual representable fine intervals; only fine states and exchanges are accepted, using the original first-order endpoint-temperature/midpoint-momentum composition and native unequal-half backward-Euler error scaling. Queries retain the actual interval midpoint, endpoint and authored continuous-source segment. The error norm includes velocity, original caloric energy and temperature changes, signed transported volume, each fluid cell's deposited viscous heat, original per-cell other heat, and separate signed/magnitude body, wall and surface work. Exact accepted source, link, ambient and reservoir overrides feed the original rate trace. Different heat mechanisms within the same cell can still cancel; no per-boundary or per-mechanism error certification is claimed.

Native boundary/contact events split intervals and execute through the original same-time `advance_to`; target events remain pending at the left limit. Any change to any native cell identity, material, mass, geometry or stored enthalpy, including nonfilm inventory, refuses the entire window. No remap, birth or advected inventory is inferred. Whole-window thermal and kinetic storage and accepted-only physical work and heat close independently. Whole-call endpoint/material ceilings include discarded work; accepted native steps also obey the original `max_substeps`. Original support and algebraic limits, cancellation and source failure remain terminal. A final step suggestion supports continuation. The film retains its original fixed thickness and density; streamwise thermal advection and evolving free-surface inventory are outside this owner.

Verification: 9 of 9 source-linked thermal-film tests passed, including five new independent adaptive regressions. The existing bounded projection uses actual film, Carreau, linear-algebra, math, Krylov and original thermal kernels with pinned current thermal inputs; it does not establish full-crate Cargo or executor/runtime integration. Targeted formatting and diff checks passed.

## Radiating spatial caloric electrothermal feedback

`advance_caloric_tensor_electrothermal_radiating` combines the original
spatial current/contact solve with the native radiating caloric backward-Euler
owner. Resistivity receives original specific energy and temperature, so
latent fraction can change both current and Joule heat at constant Kelvin.
Each source proposal retains the same opening inventory; discarded proposals
neither advance time nor deposit additional heat.

Immutable ray geometry is prepared once for the complete coupled interval
and borrowed by every thermal proposal. Original conduction, caloric inverse,
physical radiation powers and thermostat reactions remain with their existing
owners. The accepted radiation report accompanies the accepted thermal step;
its heat is already in the native boundary ledger. Applied Joule heat on
prescribed cells retains the separate equal source/thermostat-removal channel,
already included in total electrical heat.

The common electrical-source iteration retains its local source comparisons,
noncancelling joule allowance, original-law query budgets and cancellation.
The outer pass bound and native radiation/thermal bounds limit total work.
Fixed geometry, gray/diffuse/opaque surfaces, transparent gas, fixed
emissivities/surroundings and the original cell-temperature radiation sampling
remain explicit. Authored sinks must not duplicate the same radiation.

Two regression groups exercise latent/contact/thermostat balance in both
escape and radiosity modes, and nonlinear single-phase cooling against an
independent physical quartic balance. Independent scalar references and pinned
Rust formatting/parsing passed. Native Cargo test and single-threaded check
attempts were killed while compiling the original asupersync dependency, before
fs-lbm was reached. Native integration compilation and runtime execution remain
unverified; these results do not claim experimental material or radiation validation.

## Physical-stage spatial caloric electrothermal evolution

`advance_caloric_tensor_electrothermal_sdirk2` and its radiating counterpart
compose the same original current, contact and caloric heat owners at both
physical SDIRK2 stages. Separate stage samples carry electrical drives,
electrical contacts, authored thermal forcing and prescribed specific energies.
The resistivity callback receives the original stage clock, cell, occupancy,
specific energy and temperature. Both stage energy fields and refreshed local
Joule powers must converge, including isothermal latent intervals.

Every proposal calls the original two-stage heat owner from the same physical
opening inventory. The native owner alone carries first-stage heat and
prescribed-cell storage into the second row. Applied electrode energy, Joule
heat, thermostat Joule heat and local source discrepancies use the positive
weights dt*(1-gamma) and dt*gamma. Absolute local/stage defects cannot cancel;
source-combination and native tableau/reduction rounding enter the existing
joule allowance. Each stage's prescribed mask owns its thermostat channel.

The shared coupling engine supports one or two stages. It retains both applied
and refreshed electrical fields, reuses refreshed fields for the next proposal,
and moves native energy vectors without copying complete heat solves. Two-stage
current work is exactly 2*(passes+1) solves. Native thermal and radiation work
limits span both stages per proposal; the outer pass bound limits proposals,
while original source query budgets persist throughout the coupled interval.

Radiation geometry is prepared once for the whole interval. Returned stage
reports expose physical powers; their tableau-weighted integral is already in
the boundary ledger. The transient record's radiated_w remains final power.
Fixed mass/geometry, stationary DC and original gray-surface assumptions remain.
Second-order accuracy applies to smooth, resolved forcing and fixed constraints;
forcing discontinuities and sharp phase corners require interval splitting.

Three added regression groups cover stage-dependent current/contact/source
data with changing prescribed latent inventory and a late second-stage refusal,
physical radiation at both stages against independent quartic equations, and
spatial temporal refinement against an analytical exponential solution.
Independent references and pinned parsing/formatting passed. The native Cargo
attempts above did not reach these nine electrothermal regression groups.


## Grid-aligned orthotropic immersed solids in native convection

`ThermalOpenField3::set_immersed_orthotropic_thermal_voxels` and the matching
`BoussinesqOpenGrid3` entry retain the three diagonal conductivity components
of each original solid voxel. A normal solid-solid path uses the harmonic
conductance of the two corresponding half cells; a fluid interface uses that
voxel's normal half-cell and authored normal contact resistance, followed once
by the original D3Q7 reflected-pair resistance. No directional average replaces
these physical paths. Opposing faces share their normal-axis contact value.

The original immersed owner now retains each exposed link's conductance and
weight and factors the same simultaneous sparse solid/fluid-trace operator
once. Centering weighted traces on a canonical original interface preserves
uniform equilibrium and declaration-order determinism. There are no per-step
allocations, geometry searches, new inventories or separate heat advances.
`ImmersedThermalBodyState3::surface_resistance` reports the actual three paths;
`solid.resistance` remains the common scalar declaration (zero for a
nonuniform directional declaration). Capacity, source, accepted temperatures,
first-law accounts and final publication remain with the original owner.

Isotropic voxels use the same preparation and retain exact parity with an
equal-component orthotropic declaration. Original periodic parallel faces,
interior-node connectivity, node/link/factor/storage limits and transactional
installation remain active. This fixed-coefficient, cubic-grid model resolves
diagonal tensors only. It does not infer material data, rotated off-diagonal
transport, nonideal solid-solid contact, nonlinear caloric laws or geometric
surface accuracy.

Twenty original immersed-owner tests passed in the bounded scalar source
projection, including four new directional groups: independently calculated
per-face reflected heat, heterogeneous two-cell balances under all three axis
rotations and periodic parallel faces, exact isotropic parity and extreme-
capacity equilibrium, and data/budget/late-cancellation refusal with exact
retry. The exercised streaming, boundary, interface, sparse factor and thermal
transaction routines and their tests were original source. The projection
omitted the unused buoyancy module and executor wrappers, and retained only
the unused tensor carrier/shape interface with the thermal-owner error type;
no tensor was constructed. Full Cargo and public Boussinesq integration remain
unexecuted; these local balances do not establish spatial convergence or the
developing-channel Nusselt benchmark.

## Opposing rigid planes on one native liquid tick

`FreeSurface3::planar_rigid_wall_pair_candidate` and
`step_with_planar_rigid_wall_pair` admit two distinct opposing, fully flooded
domain planes on one Cartesian axis, preserving caller order and each wall's
motion and moment origin. Both original rigid-wall reflections share one
collision/streaming/VOF-conversion tick and one liquid candidate. Disjoint
plane ownership prevents a link from receiving both wall corrections.

Each original 3-by-3 affine wrench response describes the same accepted
incoming liquid. Its instantaneous cross-wall derivative is zero; subsequent
native fluid evolution carries interaction between the planes.
`PlanarRigidWallPairReport3` counts the liquid inventory once and sums actual
reflected-link impulse and mechanical work with exact wall reciprocals. Its
checked common-origin angular-impulse methods include orbital contributions,
`(wall_origin - requested_origin) cross impulse`, using double-double sums.

Existing local/relative Mach limits, endpoint flooding, mass admission and
cancellation remain authoritative. Fixed planes support tangent translation
and normal-axis spin. Normal displacement, tilt, remapping, partial wetting,
finite plate edges, solid dynamics and thermal feedback are outside this
operator. Boundary work does not certify liquid total-energy conservation.

## Source-bound opposing rigid-wall momentum in SI units

`MaterialFreeSurface3::planar_rigid_wall_pair_candidate_si` and
`step_with_planar_rigid_wall_pair_si` retain the original material/scaling and
convert the shared native tick through the existing rigid-wall SI and
pressure-reference owner. Each wall retains its raw numerical-datum wrench,
pressure perturbation load and absolute pressure load only when an absolute
reference was explicitly declared. Changing that physical reference does not
change the liquid populations or mechanical boundary work.

Aggregate impulse and work use the original physical units. Mass describes
one shared inventory; the repeated per-wall inventory observations must not
be added. Checked common-origin angular impulses include orbital moment and
exactly reciprocal wall signs. Both local wall motions and the common liquid
endpoints retain the original material scenario limits. Publication follows
both report conversions, native mass/flooding admission and the final
cancellation checkpoint. Mechanical work is not deposited as caloric heat.

Validation: all eight new native pair tests and all ten selected SI moving-wall
tests (five new pair cases and five existing cases) execute successfully against
the selected original numerical/material sources. The two older late-speed
fixtures now exceed their scenario bound from physical rest; no production
speed limits or test tolerances were relaxed. Full Cargo validation remains
blocked by the 8 GiB environment limit while compiling asupersync metadata;
unrelated async runtime registrations are excluded from the focused build.

## Two rigid wall exchanges on one thermal liquid trajectory

`MaterialFreeSurfaceEnthalpy3::planar_rigid_wall_pair_candidate_si` admits two
opposing flooded rigid wall motions through the existing source-bound thermal
owner. Authored cell heat refreshes the original temperature-dependent local
viscosity/rheology and optional Boussinesq force before one collision, both
reflections and VOF conversion. The same native mass-route trace carries
extensive enthalpy. Advected material, wetting and force-frame refresh then
rechecks both actual wall endpoints before returning the shared candidate.
The native transport entry borrows that sole candidate, without cloning the
population state or inventing a second heat-advection flux.

`planar_rigid_wall_pair_candidate_with_viscous_heating_si` additionally uses
the original held initial-state molecular stress-power sampler. Authored heat,
molecular heat, both measured wall works and body heat remain distinct. One
thermal interval advances its original accepted clock once; a discarded body
pilot can deposit no second copy of heat in the accepted state. Both wall
wrenches retain their original SI pressure reference, reciprocal impulse,
common-origin orbital/spin moment and boundary work accounts.

Constant-temperature compatibility, original zero-wall thermal replay,
heated-viscosity/buoyancy response, independent source-once molecular replay
and late source/second-wall/cancellation refusal have focused executable tests.
The first-order thermal source split retains its timestep-refinement requirement.
Fixed fully flooded planes permit tangent translation and normal-axis spin;
normal movement, tilt, changing coverage and finite edges remain unsupported.
Prepared conduction/electrical/surface-temperature owners retain their original
admission rules. Mechanical wall work is not deposited as caloric heat and
this does not claim a discrete whole-fluid total-energy balance.

Validation: four new thermal-pair tests, eight native pair regressions and five
SI pair regressions pass against the selected original numerical/source code.
A separate existing partially filled wall-field fixture refuses in ordinary
stationary transport, before any wall entry, with `native cell mass routes do
not close to roundoff`; that untouched baseline path is not claimed green.
Full Cargo/workspace validation remains outside this focused execution because
the async runtime dependency exceeds the environment's 8 GiB memory limit.

## Original phase-dependent conductivity on native caloric tensor cells

`TensorConductivityLaw::caloric_conductivity` is an optional original-energy
constitutive route. Its `TensorCaloricConductivityPoint` contains the absolute
specific energy (including its original datum), recovered temperature, storage
basis and original caloric identity. Returning `None` retains the existing
interval-integrated K(T) route and its original inverse-query behavior.
Returning `Some(TensorCaloricConductivity)` explicitly selects cellwise K(e)
on the existing HMM geometry, with a physical SPD tensor and finite symmetric
analytic dK/de. The derivative may be negative or indefinite. A Kelvin-only
entry never invents an energy or a phase fraction for this route.

Every cell/face link uses that trial cell's tensor. The original heat is
`q_f = -sum_g M_fg(K(e_c)) * (T_g - T_c)`, and the native Newton column is
`dq_f/de_c = B_f(K) * dT/de - sum_g M_fg(dK/de) * (T_g - T_c)`.
This retains phase-dependent transport even on an exactly isothermal latent
plateau. Face unknowns remain Kelvin and contacts retain the two original
half-cell resistances plus the authored contact. Signed cell derivatives are
not clamped; condensation still requires a finite positive total local
diagonal including physical storage and sinks. A nonpositive diagonal is a
bounded refusal, not a claim that every nonlinear state is solvable. No face
enthalpy, apparent heat capacity, duplicate latent inventory or new heat owner
is introduced. Backward Euler, SDIRK2, adaptive, radiating and electrothermal
caloric entries inherit this route through their shared spatial owner.

`IsotropicPhaseConductivity` borrows one original `PreparedPhaseConductivity`
without copying or resampling its paths. Constructing the view explicitly
declares isotropy. It admits strictly positive retained conductivity, checks
the original enthalpy basis and caloric identity, and preserves original
source support and analytic derivatives. It also implements caloric storage:
using the same view for storage retains restricted bounds and the source's
selected one-sided dT/dh at an endpoint. Original density and entropy remain
with the original caloric curve. Kelvin-only queries are refused because
temperature alone cannot select latent inventory. The hot path uses the
prepared owner's bounded interval lookup, without table allocation or source
selection. This cellwise K(h) model does not claim an exact Kelvin Kirchhoff
integral or replace the separate material P1 transport owner. Geometry and
reference mass remain fixed; no expansion or moving geometric fronts are
inferred.

All 19 tests in the two actual numerical leaf modules pass with the pinned
Rust compiler. They include full-tensor analytic/finite-difference Jacobians,
latent and sensible derivatives, zero-dK/de parity, and signed Schur recovery
against a complete independent 7-by-7 solve. Four new integration groups
cover constant-Kelvin phase heating and datum shifts, BE refinement, spatial
contact transfer with independent BE/SDIRK equations, source/cancellation
refusals, and an actual MaterialCard/prepared-source solve crossing a latent
conductivity knot. Their independent scalar references and pinned parser
checks pass. Native production and integration-target type checks also pass;
the integration-execution resource limit is recorded below. Numerical-leaf
execution is separate from full native integration evidence.

### Thermal feedback in resolved film mechanics

`ThermalCarreauFilm::advance_with_thermal_drive` and `advance_adaptive_with_thermal_drive` resolve body acceleration, wall velocity and surface traction against every actual original caloric endpoint candidate, including derivative probes and the independent final evaluation. Each proposal restarts the same accepted momentum history. Heat/interface resolution receives the exact candidate drive and its owned source payload. The accepted report retains that sample, its interval query, explicit thermal endpoint time and unchanged payload; the adaptive report retains only accepted fine-step samples in chronological order. Payloads need not implement Clone.

The prescribed `advance` path retains its held drive, and prescribed `advance_adaptive` still samples its source exactly once per native interval at the authored midpoint. Thermal derivative probes do not multiply prescribed callback queries. All modes share the original nonlinear residual and adaptive controller. A separate drive-query ceiling and actual count span the complete operation, including discarded comparisons and failed queries; existing endpoint, material and original native-step ceilings remain independent. Callbacks must remain repeatable and free of physical side effects.

Original body, wall and surface work remain separate signed mechanical exchanges. Only the original resolved viscous loss enters caloric storage; source payloads cannot replace physical accounting. Endpoint-held feedback retains the existing first-order coupled composition, event/left-limit semantics, local state and transfer error indicators, original source/rate/support limits and independent energy gates. A source, budget, physical or cancellation refusal publishes no partial window. No actuator model, thermal expansion, changing mass or new constitutive fit is inferred.

Four G1/G3/G4 regression groups are authored for an independent one-cell momentum/heat root with the original positive Arrhenius viscosity, transient refinement against continuous equations across a scheduled reversal, exact prescribed/feedback parity and callback counts, and late query-budget/source/cancellation rollback with retained source identity. Independent scalar reference calculations and pinned Rust 2024 parsing/format checks passed. The new Rust tests remain uncompiled and unexecuted because the prior bounded material/solver build setup was lost with the execution environment; no full-crate Cargo or runtime result is claimed.
## Material enthalpy in original immersed-body convection

`ThermalOpenField3::set_immersed_body_caloric` and the matching
`BoussinesqOpenGrid3` method bind one newly installed independent immersed
body to its original immutable `EquilibriumEnthalpyPhaseCurve`. Explicit
fixed mass, enthalpy span, temperature origin/scale and fluid-cell capacity
convert the SAME native thermal coordinate to SI. Binding must reproduce
the installed initial temperature exactly; rebinding, post-exchange binding
and any conducting body/voxel graph are refused before mutation.
All caloric bodies bound to one field must share identical temperature
origin/scale and reference fluid-cell capacity, preserving one physical
interpretation of the shared native fluid traces and summed heat.

The original body now retains one accepted `ImmersedCaloricState3`, including
the actual material/phase identity. Its native temperature is a derived cache;
the initial fixed capacity is inactive declaration data when this state is
present. The retained Arc preserves the original curve and source receipts.
The original material query supplies exact affine-Cp sensible calorics and
isothermal latent plateaus without fitting, extrapolation, phase smoothing or
an artificial apparent heat capacity. Fixed mass never follows chart density.

The same post-collision interface traces enter one bounded backward-Euler
balance in specific enthalpy. Safeguarded Newton/bisection uses the original
T(h) derivative, including exactly zero on a latent plateau. Only the required
side of the admitted source bracket is queried; a bulk-fluid temperature
outside the solid material span does not demand an irrelevant material query.
Actual reflected incoming-minus-outgoing heat and actual representable
m*(h_new-h_old) independently pass the original local and joint fluid/solid
balance checks, as well as the authored physical caloric tolerance. Source
heat is counted once; no latent-heat source is added. Cancellation is polled
at each bounded iteration/source bracket and the original streaming work.
Native populations, clocks, caloric state and receipts publish together only
after success. There is no per-step allocation or new thermal trajectory.

This capability is for isothermal retained bodies with fixed geometry and
fixed authored interface resistance, including a lumped phase-change store.
It does not remove the obstacle when liquid fraction grows, model liquid
motion, resolve internal gradients, update k(T), infer volume from density,
or close expansion/pressure work. Those approximations require caller
justification. Spatial nonlinear caloric voxel graphs and the developing
channel Nusselt benchmark remain outside this delivery.

Twenty-eight original immersed-owner tests passed in a bounded scalar source
projection, including eight new groups: constant-capacity parity and actual
original affine-Cp calorics against an independent quadratic, exact melting
and freezing plateaus and phase crossings, SI/native energy closure,
source-domain/work/representability refusal, and cancellation with exact retry.
The projection compiled complete original fs-blake3, fs-obs, fs-math,
fs-sparse, fs-qty, fs-evidence, fs-ivl and fs-matdb dependencies. The original
material phase core, source ingress, source-receipt types, caloric derivatives
and identity methods were retained, omitting only unrelated phase children.
The exercised native scalar kernels, immersed owner, sparse paths and tests
were original. The unused native buoyancy wrapper and tensor constructors
were omitted; the retained unused tensor carrier/shape methods used the
thermal-owner error alias and no tensor was constructed. Full Cargo and the
public Boussinesq wrapper remain uncompiled in this environment. These tests
establish local implementation behavior, not spatial convergence, empirical
PCM validation or a complete multiphysics pressure-work account.

## Evolving porous laws through native open flow and heat

`ResistanceOpenGrid3::{advance_serial, advance_pooled}` supplies authored
next-endpoint resistance, material axes, skeleton velocity and body force.
Collision retains the accepted force and drag; native streaming, complete
reservoir traces and their interior donors resolve the next law at their own
physical velocity. Mixed wall/open seams still require exactly zero drag.
The held `step_*` paths retain their original arithmetic and meanings.

`ResistanceOpenEndpointExchange3` combines the original shared
`ResistanceEndpointExchange3` with actual native reservoir transport. Wall
links, wall work and periodic angular transport retain their existing
measurements. Endpoint-power dissipation is nonnegative; impulse work and
its signed temporal quadrature remainder remain distinct. A finite explicit
budget bounds the sum of absolute per-cell relative work remainders before
publication, so spatial cancellation cannot hide an unresolved source history.

`advance_serial_with_thermal` and its tensor-diffusion counterpart reuse the
original D3Q7 heat transaction. Opening advection/source history uses accepted
physical velocity; the closing source half-kick uses the new law's resolved
endpoint velocity. Laws, populations, forces, drag, thermal walls/history,
both receipts and clocks publish together. Failure leaves the accepted fields
retryable. `last_exchange` and `last_endpoint_exchange` identify the latest
accepted held or evolving interval; the other selection is None.

No new collision, porous root solver, thermal transport kernel or per-cell
scratch allocation is introduced. The caller supplies a complete next law
field; this does not infer a temperature-dependent permeability closure,
drag-heat partition, porosity evolution, skeleton inertia or transported
total-fluid-energy balance. Fixed-geometry and low-Mach limits remain.

Seven focused G0/G1/G3/G4/G5 regression groups cover independently manufactured
changing tensor/skeleton/Forchheimer endpoints and common boosts, complete
reservoir/wall/voxel/periodic balance, serial/pooled equality, held-law bit
parity, invalid seams/laws and temporal budgets, actual thermal force epochs,
finite-wall heat balance and late cancellation/refusal with deterministic retry.

## Proposed boundary targets in native flow/heat publication

Both `OpenBoundaryGrid3` and `ForcedOpenBoundaryGrid3` provide
`step_serial_with_thermal_and_open_boundaries`. The existing boundary admission
prepares targets without publishing them; the original collision, streaming,
force-aware reconstruction and D3Q7 heat operators use one shared final gate.
`OpenThermalDiffusion3` selects the original borrowed scalar or tensor law.
The forced entry holds its current field or swaps an explicit endpoint buffer
only after acceptance. Existing thermal entrypoints use the same implementation.

Proposed wall targets participate in thermal ownership admission, including the
moving-finite-storage refusal. Retained profiles, finite wall/conduction history,
heat flux history, actual opening/endpoint velocities and separated momentum/heat
accounts keep their original semantics. Scratch may change after refusal; accepted
controls, physical fields, sources and counters do not. No heat is inferred from
body force or mechanical wall work, and the original temporal-order limits remain.

A focused native regression stages real momentum before closing-heat/final-cancel
refusal, checks unchanged finite walls and force buffers, then reproduces an exact
retry. Coupled tests exercise both lattices and diffusion families. Pinned parsing
and scoped diff checks passed; these tests were not compiled or executed here.

## Native conversion roundoff in extensive thermal transport

The mass-route trace carries the measured residual of each conservative VOF
conversion only after the original donor-density-scaled admission succeeds.
The absolute amount is attached to its converting source node and accumulated
locally. Thermal mass closure adds that original allowance to its existing
throughput-roundoff bound. An almost empty source therefore does not lose the
roundoff scale of the finite donor inventories used to repay its conversion
debt. No global mass allowance, unit floor, fictitious route or mass repair is
introduced.

Native mass transfers, momentum transport, source laws, matrix coefficients
and solver arithmetic are unchanged. The shifted enthalpy equation retains
its explicit `(supplied - removed) * reference` contribution, so admitting an
existing mass residual does not silently remove its energy consequence.
Unaccounted defects and wrong-node allowances continue to refuse.

Focused original-source execution passes the existing four-tick partially
filled stationary/moving-wall thermal regression unchanged. Two new tests
check the reference-energy term at mass scales 2^-60, 1 and 2^60, local absolute
accumulation, missing/wrong/insufficient allowances, larger defects and
nonfinite/overflow refusal. This is focused native numerical validation, not a
full Cargo feature-graph or repository-level gate.


## Adaptive scheduled caloric electrothermal evolution

`advance_caloric_tensor_electrothermal_adaptive` advances a complete original
energy/current window to an exact target time, with optional physical gray
radiation. `CaloricElectrothermalForcing` retains each stage's thermal forcing,
electric setup and electrical contacts. Forcing and resistivity receive the
same original `TensorCaloricForcingSample`, including an authored transition's
left-limit marker. The fixed prescribed mask, event clipping, stage clock,
embedded controller and compensated heat ledger remain with the existing
native adaptive owner. Each coupling proposal restarts the same physical
opening energy; rejected temporal/coupling proposals publish no physical heat.
Original source objects and their query budgets persist across the whole run.

Temporal admission checks both original caloric increments and each cell's
applied Joule-energy discrepancy. The latter compares the actual SDIRK2
weighted stage heat with its first-order embedded estimate using the original
scale-safe error owner. Prescribed cells participate: a thermostatted conductor
with constant energy still refines time when its electrical heat changes.
Electrical source error cannot cancel between cells. Relative caloric scaling
uses changes, preserving the original energy datum. These local asymptotic
estimates do not certify global error or every independent external heat
channel. Declare source jumps; sharp material corners can reduce time order.

Radiative geometry is prepared once for the whole window. The result retains
small accepted coupling/radiation receipts and only the last accepted
interval's electrical fields. Accepted electrical work, Joule heat, thermostat
heat and noncancelling coupling defects accumulate in double-double accounts.
Native thermal closure and the accepted coupling-defect sum each face the
whole-window energy allowance. The thermostat electrical receipt is a separate
source/removal account; it is not deposited again in native caloric storage.

Completed-trial work includes every successful thermal proposal inside each
completed coupling trial, including trials rejected for time accuracy. Incomplete
failed trials remain separately bounded by attempts, coupling passes and native
solver ceilings; their partial work is not claimed as completed work. Late
source, budget, energy or cancellation refusal discards the whole private prefix.
This remains fixed-inventory, fixed-geometry stationary DC with the original
gray-surface radiation model and immutable original material laws.

Four integration groups cover latent exponential feedback and datum invariance,
changing heat through a fully thermostatted conductor, scheduled current and
resistivity transitions with late refusal/cancellation, and radiating cooling
against an independent nonlinear ODE. The existing electrothermal sensible
fixture now uses the actual single-phase constructor and explicitly typed
conductivity arrays. All 24 tests in the three actual numerical leaf modules
pass with the pinned Rust compiler. The native production library and all four
named integration targets (`tensor_caloric_phase_conductivity`,
`tensor_caloric_electrothermal_adaptive`, `tensor_caloric_electrothermal`, and
`tensor_caloric_radiation`) pass `cargo check --locked` with enthalpy coupling,
authentic locked sibling sources, one compiler thread, no debug information and
MIR optimization disabled. Actual execution of those four targets was blocked
by SIGKILL during asupersync code generation under the 8 GiB memory limit,
before any test binary ran. Native type checking and the 24 numerical-leaf
passes are not a claim that those integration tests executed.
