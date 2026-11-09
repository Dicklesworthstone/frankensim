<!-- fs-couple contract, part 5 of 8; a higher number is newer. Index: ../CONTRACT.md. -->
## Finite solvent heat during gas absorption and exsolution

`cavitation::dissolved_gas::advance_thermal_mixture` closes the original
finite-inventory Henry/two-film kinetics with one finite incompressible solvent
and a common gas/liquid temperature. The initial-temperature mass substep is
the original simultaneous owner, including opposite species transfers. One
existing safeguarded scalar solve then balances its actual internal-energy
change with original gas, dissolved-partial and pure-solvent sensible increments
and the prescribed external heat. Zero external heat makes the combined pair
adiabatic. The isothermal substep's thermostat accounts are an explicitly
reported numerical decomposition, not additional physical heat reservoirs.

All caloric differences use the original source-integrated increments; gas
formation datums are not subtracted as large absolute storage totals. Gas mass,
composition and volume survive the temperature proposal exactly, with pressure,
transport and acoustics reconstructed once by the original material owner.
The same Henry/diffusivity claims continue at the accepted T/p/concentration,
including actual interfacial concentrations, shared mixture dilution and total
interface-pressure admission. No nested gas-temperature inversion occurs in
the coupled heat root. Failed source, energy, work-budget or final cancellation
checks publish none of the gas, liquid clocks or refreshed source laws.

This is a first-order conservative split with explicit time and temperature
step restrictions; energy closure is not a kinetic trajectory-error bound.
`dissolved_gas::trajectory::advance_thermal_to` resolves a held-power window
by comparing one original full step with two original half steps, accepting
only the latter. Original H(T), liquid and gas diffusivities, pressure and
calorics refresh at every native endpoint. Componentwise local errors cover
temperature, pressure and every gas/dissolved species; relative scales use
changes, so a large background cannot hide a small evolving coordinate.
Actual representable half durations scale the first-order estimate, including
unequal halves at large clocks. This remains an asymptotic local estimate,
not a global accuracy or empirical uncertainty bound.

Only time-error and typed temperature-increment failures permit shorter
retries; original source/regime/caloric refusals retain their authority.
Native evaluations are explicitly bounded. Accepted physical increments alone
enter compensated whole-window heat accounts; original/final species
inventories independently enforce fixed per-species and total mass budgets.
These conservation tolerances are never multiplied by the number of steps.
Gross absorption and exsolution retain transfer reversals in bounded memory,
alongside the final native receipt. Cancellation or failure discards the entire
window, including accepted internal prefixes. Five additional G1/G3/G4 groups
cover tolerance tightening against an independent nonlinear RK4 trajectory,
warming-induced exsolution with original gas-film sources, sealed finite heat,
opposite-species ordering, unequal representable halves, and transactional
minimum-step/work/source/final-cancellation refusals.

Five G1/G3/G4 regression groups cover independent constant-Cp absorption and
exsolution endpoints, the sealed-mass native heat limit, thermal source refresh,
16/32/64-step refinement against an independent energy-constrained RK4 oracle,
opposing-species canonical results, and late source/cancellation rollback.
The native tests have not been executed in this session because the pinned
runtime dependency build was killed by the shared memory limit.

The solvent chart must remain liquid at constant density and solute partial
volume is negligible. Independent Henry, partial-caloric and gas-film sources
must apply over the complete pressure/composition/temperature interval and be
physically compatible; this numerical composition does not certify that data
compatibility. Finite interfacial thermal resistance, solvent evaporation,
moving bubbles, spatial transport, reactions and entropy-production claims
remain outside the common-temperature dilute closure.

## Joint sourced dry friction and finite contact heat

A `DryFrictionBinding::new_temperature_source` retains the original
`SourceThermalCoulomb` claims and ordered material pair. Each heated evolving
contact contributes an interface-temperature unknown to the existing finite-cell
enthalpy solve. At every trial endpoint the same pinned coefficient claims are
queried at that candidate temperature, then the original simultaneous mechanical
solve supplies traction, slip, normal pressure and dissipated friction work.

For each massless contact plane the additional residual is
`sum_i(dt * G_i * (T_contact - T_cell_i)) - D_friction = 0`.
Signed fluxes enter the original cell enthalpy account. The temperature coordinate
uses `dt * sum_i(G_i)` only as a numerical unit scale, not an added heat capacity.
Interface, cell and total energy residuals must pass before coefficients, source
receipts, mechanics and finite thermal histories publish together. Failed trials
restore the accepted laws and load-line caches; release disconnects conduction
and friction heat without querying a positive-pressure source at zero load.
Constant laws retain their original analytic interface-temperature path.

This route admits temperature-dependent Coulomb coefficient values. Actual
normal pressure and incoming, average and endpoint slip remain applicability
checks; pressure-, speed- or additional-history-dependent coefficient responses
refuse. The native material-string input retains this source through ribbon area
partition and requires a finite dry heat route for every evolving physical
contact. No recording fit, invented reservoir, missing rosin data or transported
third-body layer is supplied by this coupling.

## Sourced directional force/contact/drift for mixed articulated mechanisms

`friction_heat::directional::articulated::step_source_mixed_articulated_elliptic_contacts`
connects original material-interface histories and coefficient laws to the native
loaded mixed-tree step. Every geometrically active source is sampled once from
actual post-load-kick tangent kinematics; inactive sources remain unqueried.
Original ordering, claim pins, coefficient role, texture and accepted travel
coordinates retain the existing source sampler. Native simultaneous directional
contact, free/prescribed drift, common clocks and complete-path admission remain
the original mechanical owners.

One candidate carries actual endpoint trees and source receipts only after all
native source, drift, endpoint, momentum, work and final-cancellation gates.
The caller adopts them together. Complete native boundary work includes driven
supports as well as moving planes; load support work, drift support work and
numerical force/contact splitting remain separate. This first-order held-source
composition does not integrate source travel or infer pressure, wear, thermal
partition or heat from contact/splitting loss. Conditions, support-path and
geometry callbacks must stage their external state until complete acceptance.

Focused original-source regressions specify analytic post-kick speed, both joint
recoils and actual travel, driven support work, common boosts, inactive histories,
typed source/path refusal, terminal cancellation and deterministic retry.
The source and tests received scoped Rustfmt parsing and API review; compilation
and execution remain unverified in the constrained editing environment.

## Sourced directional friction in articulated impacts

`friction_heat::directional::articulated::resolve_source_articulated_elliptic_contacts`
connects the original ordered material-interface claims to simultaneous contact
response of free-floating jointed bodies. It reuses the existing source binding
and preimpact coefficient sampler, with the native link-point kinematics,
joint constraints, base recoil and full cross-contact mobility.

`resolve_source_mixed_articulated_elliptic_contacts` applies that same source
sampler to free and prescribed-base trees through the native mixed callback
owner. Its generic candidate retains the actual native response type, including
full support reactions and signed work. Total `boundary_work` includes both
plane and base-support work; under nonunique impulses only the total is
identified by the native work witness. The original free API keeps its default
response type. Mixed regressions specify both joint recoils, driven support
work, pinned source histories and late angular-momentum refusal with retry.
Additional native source tests cover full-tensor anisotropy, rotated material
frames and exact free-slip directions. These new tests have not been compiled
or executed in the constrained editing environment.

Each geometrically active contact is sampled once in canonical contact-ID
order using its actual incoming tangent speed and caller-owned local
conditions. The exact normalized frame, both physical witness positions and
original input index preserve source/history mapping; inactive histories stay
untouched. Original claim pins, accepted path coordinates, texture and complete
static/kinetic continuation state remain attached to the sampled candidate.

One unpublished result contains every native mechanical response and the
complete active source candidates. The adapter borrows accepted inputs and
avoids cloning complete source histories into a second staging array. The
surrounding event owner must publish all returned body responses and source
histories together. Source refusal, cancellation and late mechanical admission
failure return no partial candidate; callbacks must stage their external state.

Four G1/G3/G4 regression groups specify independent prismatic-joint/base-recoil
and energy references, heterogeneous sources and caller permutation, inactive
histories, later-source refusal/cancellation and late momentum refusal with
deterministic retry. These synthetic references do not establish experimental
calibration.

This is a held-source instantaneous impact law with an explicit coefficient
role. It does not infer pressure from impulse, advance finite travel, allocate
impact loss to heat/wear or automatically transport texture between changing
frames. The caller supplies local conditions and uses the existing objective
history transport when the contact frame changes.

## Sourced simultaneous directional rigid impacts

`friction_heat::directional::impact::apply_source_rigid_elliptic_impacts`
binds each original ordered material interface to the existing simultaneous
rigid-body contact owner. Every geometrically active source is sampled once at
its actual incoming tangent speed and caller-owned local thermal/pressure
conditions, before normal and tangent impulses are solved together. Original
body masses, rotated inertia, material arms, normal/tangent coupling and
shared-body cross-contact response remain in the native mechanical assembly.

Source callbacks receive the exact normalized contact frame, both actual
witness positions, stable contact ID and original input index. A prescribed
plane's source point is the normal projection of the body witness, not an
arbitrary plane anchor. Canonical contact-ID ordering retains the original
binding/history mapping; inactive contacts neither query nor advance sources.

The source sampler retains claim pins, accepted accumulated distances, texture,
native slip and complete static/kinetic continuation state. Source receipts,
body momenta and candidate histories belong to one accepted event. Source
failure, cancellation, native non-admission or a late physical momentum/work
audit leaves both public state arrays unchanged. Callbacks must stage any
external state themselves. A changed tangent frame requires the existing
objective history transport before the next event.

Four new G1/G3/G4 groups cover analytical off-centre normal/tangent coupling,
heterogeneous three-body friction with independent impulse/endpoint references,
caller permutation and inactive histories, later-source refusal/cancellation,
and late angular-balance refusal with deterministic retry. All eleven sourced
sphere and general-impact groups pass. The stronger-anisotropy shared-contact
fixture retains explicit unresolved-determinacy refusal and complete rollback;
positive coefficients alone do not imply native uniqueness admission.

This is a source-held instantaneous impact law, with an explicit static or
kinetic coefficient role. It does not infer pressure from impulse, a fictitious
duration, finite travel, collision discovery or thermal allocation. Impact loss
remains unallocated mechanical loss. These synthetic source references do not
establish experimental calibration.

## Energy/momentum stepping of nonlinear finite bodies

`FiniteSolidBody::step_energy_momentum_candidate` opts into the prepared
solid's objective discrete-gradient force and analytic tangent while retaining
the native midpoint/generalized-alpha Newton-FGMRES solve, source material,
consistent reference mass and fixed supports. The operator reconstructs the
endpoint only inside midpoint force/tangent callbacks; actual-state validation
keeps its original meaning. Endpoint differentiation supplies the required
factor of two in both tangent action and preconditioning.

Held nodal and body dead loads retain displacement-conjugate work. Support
reactions use the same algorithmic force and full consistent-mass acceleration.
Pressure-bound bodies use the explicit conservative-pressure energy/momentum
entry for each held interval. The ordinary midpoint path remains available. Stable midpoint
velocity times the mass-applied velocity increment evaluates kinetic change;
the measured energy gate retains potential-work and nonlinear error instead
of converting it to heat. Material failure, failed convergence, energy refusal
or cancellation publishes no candidate state.

Four production stress/tangent groups and two native 250-step, h=0.02 free P1
tetrahedron trajectories were executed. For the Mooney--Rivlin fixture, maximum
energy drift changed from 3.9257e-3 J with ordinary midpoint to 1.0419e-13 J;
the original source-card-compatible Ogden alpha=2 fixture changed from
3.7880e-3 J to 5.2625e-14 J. Linear and angular momentum drifts remained at most
3.1e-14 in these fixtures. They are numerical checks, not measured material
validation or general performance results.

Public integration regressions cover free tumbling, supported dead-load work
and atomic stalled/cancelled refusal. The current focused native qualification
is stated with conservative spatial pressure below; full workspace Cargo
compilation remains blocked by unavailable sibling repositories. No runtime
stubs or substitute numerical solver were introduced. The bounded
midpoint method still requires temporal/mesh refinement and does not supply
frequency-unwarped modal flow, nonconservative per-face pressure work or
completion of the broader instrument-realism program.

### Uniform pressure with energy and angular-momentum balance

For explicit signed uniform pressure,
`FiniteSolidBody::step_energy_momentum_candidate_with_pressure` retains the
original pressure-bound body and consumes the prepared discrete gradient of
`W_material(F)+p det(F)` and its analytic tangent. Force callbacks, tangent
actions and support reactions use the same combined algorithmic force; ordinary
face-force subtraction is bypassed on this path, so pressure is applied once.
The original pressure surface owns closure/reference admission and independent
work `-p (V_after-V_before)`. Material storage and pressure-source work remain
separate. Zero pressure follows the unpressurized arithmetic exactly; the
ordinary midpoint per-face pressure API retains its existing behavior.

The conservative algorithmic volume load is consistent under time refinement,
rather than an exact finite-step boundary-traction history. Uniform pressure
has no resultant force or torque on the free closed body; its algorithmic
contribution preserves the same midpoint angular balance while doing the
correct volume work. Pressure may change between intervals; no fluid storage
or pressure-dynamics model is inferred.

Six native stress/tangent groups pass. A 125-step free-tumbling fixture at
p=3 and h=0.02 has maximum K+U+pV drift 7.11e-15 J, linear-momentum drift
5.55e-16 and angular-momentum drift 2.24e-15. Across 470 native steps, the
resolved h=.01/.005 comparison against h=.00125 over .2 seconds reduces
state error by 4.1395. The earlier coarse .04/.02 comparison was outside the
resolved asymptotic regime; energy closure alone does not certify trajectory
accuracy. These are numerical fixtures, not material calibration.

Three public regressions cover pressure-driven tumbling/refinement, changing
signed pressure against an independent breathing/work/reaction oracle, and
zero-pressure/atomic refusal. These original Cx-bound regressions now pass
in the focused native closure described below; full workspace Cargo
qualification remains blocked by missing sibling repositories.

### Conservative spatial pressure in native energy/momentum motion

`FiniteSolidBody::step_energy_momentum_candidate_with_pressure` also accepts
`TetPressure::Affine` and `TetPressure::FreeSurface`. Existing scalar callers
retain uniform-pressure behavior. The prepared native product discrete gradient
supplies force, endpoint tangent and reactions to the same generalized-alpha /
Newton-FGMRES owner; the original surface supplies actual endpoint pressure
work, including wet/dry transitions. Material storage, consistent reference
mass, source admission, fixed supports and transactional refusal retain their
existing owners. The pressure is applied once, and no numerical defect becomes
heat or inferred fluid storage. `PerFace` remains outside this conservative
entry, including an all-zero face field.

The prescribed spatial field is held over each interval. Buoyancy supplies
momentum along its gradient axis; when otherwise unforced, transverse linear
momentum and angular momentum about that axis are conserved. This does not
advance the surrounding liquid or its free surface. Finite-step loads approach
native surface traction as the interval shrinks; energy balance alone does not
certify trajectory accuracy.

The complete original source-bound material/Cx/prepared/time/coupling path was
compiled and executed in a focused native-source closure. Two new regressions
pass: 30 actual Newton intervals for each spatial field, at least two wet/dry
vertex crossings, independent total potential and pressure work, admissible
momenta, receipts and unchanged accepted state after refusal/cancellation.
All six existing energy/momentum regressions also pass, including uniform
pressure refinement, support reactions and zero-pressure identity. No material,
mass, cancellation or integrator substitute was used; full workspace Cargo
qualification remains blocked by missing sibling repositories.

## Material-bound finite bodies under spatial follower pressure

`finite_solid::FiniteSolidBody` accepts the original pressure owner's
`TetPressure::PerFace` through its existing pressure constructor and candidate
step. One held differential pressure per canonical exterior face permits
localized patches and separate pressure regions. The original material laws,
reference mass, supports, initialization and implicit-midpoint time owner are
unchanged; the Newton/FGMRES operator retains the generally nonsymmetric pressure
tangent without symmetrization. Every interval still requires an explicit load.

Uniform pressure retains endpoint-volume work. Face fields use the pressure
owner's exact straight-displacement-segment work, including rigid translation
of loaded patches. This is a path-dependent approximation to a curved dynamic
trajectory and requires temporal refinement. Numerical integration defects
remain signed energy-admission checks, never physical heat. Invalid pressure
coverage, material/geometry refusal, failed energy admission and cancellation
preserve the accepted body. G1 material-backed axial-load tests use independent
acceleration/midpoint/work formulas for two material cards; G3 checks uniform
and zero-load compatibility; G4 checks failure and deterministic retry.


## Source-owned moving-sphere collision and support

`friction_heat::directional::sphere_plane::advance_source_sphere_plane` binds
original ordered material sources to the existing sphere event engine. Support
uses source-consistent actual travel and native segment replay; each impact uses
an explicit static or kinetic source coefficient ellipse sampled at its actual
preimpact speed. Collision does not advance accumulated slip/path. Full source
static/kinetic coefficients remain available to subsequent sustained contact.

Caller-owned temperature/pressure and complete queries refresh at real event and
support boundaries. Each support bracket is supplied for that actual interval;
no pressure, area, interface temperature or collision duration is inferred.
Body response, source history and bounded chronological impact receipts publish
together only after the original motion/work gates and final cancellation check.
Later support refreshes do not overwrite retained impact claim receipts.

Sustained traction-slip work and total unallocated impact loss stay separate. The
plane owns the texture frame; this point-contact adapter supplies neither a
rotating finite thermal contact map nor an impact-to-heat allocation. Independent
source, native sliding and attitude refinement retain their existing meanings.
The source-sphere impact tests cover an analytical speed-dependent impulse and
spin, impact followed by finite arrest and static hold, late source/cancellation
rollback, and a refused collision after ballistic flight. These are authored
synthetic reference cases, not experimental validation.

## Accelerated moving-contact friction heat

`friction_heat::moving::trajectory::QuadraticContactPath` supplies actual
material-grid and footprint placements at the start, midpoint and end of a
thermal interval. `MovingFrictionHeatStep::new_quadratic` and the linear
constructor share the original finite-cell thermal, massless-interface,
ambient-coverage, entropy and acceptance owners.

Every quadratic edge crossing splits the swept overlap into polynomial
intervals. Area and its first normalized-time moment are integrated with
bounded work per cell pair, retaining the original sweep and patch budgets.
Whole-path containment includes extrema between supplied samples, so agreeing
endpoints cannot conceal an excursion beyond the material surface.

`affine_power_work_j` consumes nonnegative physical friction powers at the
interval endpoints and distributes their integral by actual spatial
residence. Held traction and affine slip velocity supply this power history.
For a unit footprint moving as `x=t²` with `P=2t`, mean patch areas are
`[2/3,1/3]` but deposited work is `[1/2,1/2]`. Mean area still owns thermal
conductance and ambient coverage; it is not substituted for power-weighted
work. Area, first-moment and total-work checks admit only roundoff correction.

Native traction-to-finite-storage tests use the original friction, enthalpy
and entropy paths, with analytical cell energies and failure/retry checks.
Geometry tests cover acceleration, reversal, unsampled extrema, linear
reduction and common accelerating translation. Exact residence integration
does not remove endpoint thermal-coefficient time error. The route remains
a fixed-width planar strip with uniform spatial friction power; curved or
rotating spherical material surfaces require their own admitted geometry.

## Driven directional friction on the original finite thermal timeline

`friction_heat::directional::driven::apply_source_directional_driven_heat`
couples explicit external forces/torques, the full-tensor original rigid pair,
source-owned directional static/kinetic friction and the existing finite
enthalpy cells. The driven and unforced entry points share one material,
geometry, source, nonlinear and heat transaction. They do not duplicate a
thermal store, material law or interface heat-allocation rule.

Original principal claims are evaluated at the solved surface/interface
temperatures, load over actual contact area, mean unsigned speed and each
owned accumulated path. Resting/Sliding mode is accepted state alongside body
momenta and source history. Temperature-dependent static capacity can trigger
breakaway. Finite-cell heat receives only physical frictional work; supplied
mechanical work remains separate even when a stuck contact allows common-body
acceleration. The joint balance is
`delta(K + H_cells) - Q_external - W_mechanical`.

`advance_source_directional_driven_heat` reuses the original adaptive
full-step/two-half-step controller and retains original accepted transport
payloads. Converted distance-search bounds are rounded inward when needed so
reconstructed speeds stay inside the admitted closed source domain. The whole
requested span publishes bodies, source history, mode and thermal state once.
Compensated drive-work totals enter the independent
whole-span mechanical and combined energy gates. A late source, numerical,
budget or cancellation failure leaves every original owner unchanged.

Phase tests cover static common-body work without heat, independently
calculated breakaway work and body/cell energy, original temperature/speed/path
feedback, zero-load acceleration with conduction, and atomic refusal/retry.
Driven march tests cover analytic acceleration and arrest/hold, refinement of
a changing source law, and rollback after successful private phases.

Poses, inertia, arms, frame, load and explicit wrenches remain held during a
span; normal recoil is reported, not silently constrained. General curved
native motion is a bounded first-order discretization with its own
`maximum_step_s`. Outer phase refinement controls source/thermal sampling and
can share a fixed inner discretization error; both axes require refinement.
These local estimates do not certify global accuracy, discontinuous branch
uniqueness, spatial heat resolution or material calibration.

## Adaptive time refinement of sourced directional friction and heat

`friction_heat::directional::trajectory::march::advance_source_directional_heat`
advances a held-geometry span using the original coupled phase. Every comparison
starts one full step and two half steps from the same accepted checkpoint. It
compares original linear/angular momenta, each finite cell's energy, unsigned
path, net slip and physical friction work in explicit physical error scales.
Relative scales use increments from that checkpoint, preserving independence
from background momentum and thermal reference energy. Only the actual two-half
result is accepted; no extrapolated momentum, heat or material state is created.
Time refinement updates the original source coefficients and thermal endpoints
as bodies slow down, stop and continue exchanging heat.

An explicit mean-speed domain is intersected with the supplied phase-distance
domain using each represented duration. Smaller timesteps cannot accidentally
turn a fixed trial distance into an out-of-domain source speed, and no bound
based on initial Euclidean speed is assumed. Invalid/nonfinite domains refuse.
Source/path/thermal/geometry refusals propagate rather than being reclassified
as successful refinement. Native force, quadrature, source and Newton limits
remain active; separate comparison and accepted-phase budgets bound total work
and retained accepted payloads. Cell comparisons poll cancellation in bounded
chunks. Source callbacks evaluate private candidates and must not mutate
external physical state.

The complete requested span is atomic, including a failure after private
successful steps. Only the chosen two bodies, original source history and
original thermal timeline publish at the end. Compensated physical-transfer
totals and direct origin-to-end kinetic/cell-energy changes check whole-span
mechanical, thermal and combined closure against an absolute energy tolerance.
Every owned accumulated-path coordinate is independently compared with total
curved travel, limiting accumulated source-root and rounding defects. Accepted
phase reports retain their exact start/end times and source transport payloads.

This is local step-doubling error estimation for a fixed-pose, fixed-load,
endpoint thermal/source discretization. It does not certify global error,
source-root uniqueness, spatial thermal convergence or discontinuous-law
convergence. This unforced entry point adds no drive, normal dynamics or thermal
event scheduler. G1/G3/G4 tests include three time rungs against an independent continuous
temperature-dependent friction ODE, adaptive arrest followed by conduction,
source-speed domain admission and whole-span rollback/retry after late source,
cancellation and work-budget failures.

## Sourced directional arrest coupled to finite heat

`friction_heat::directional::trajectory::apply_source_directional_trajectory_heat`
joins the original source-owned continuous directional relaxation, full-tensor
rigid-body impulse owner and existing `FrictionHeatStep` residual. The native
trajectory provides integrated impulse, curved net displacement, unsigned path,
finite-arrest time and physical traction-slip work as separate quantities.
Average force times net displacement is not substituted for curved work.

Each thermal trial starts from the same accepted source and two body states.
Original pinned principal claims are queried at the solved finite surface
temperatures, optional common interface temperature, load/actual-area pressure,
phase-average unsigned speed and each accepted source path plus actual travel.
The same original caloric objects and located contact face are required. The
common zero-capacity interface derives signed heat into the existing finite
cells from their geometry/conductivity; no heat fraction or parallel store is
created. Physical friction enters once, while hot-to-cold exchange continues
after sliding has stopped.

Exactly the selected two bodies, source history and original thermal timeline
publish together after native impulse/momentum/work checks, finite-cell/interface
residuals, independent Delta(K+H)-Q_external admission and cancellation. Source,
domain, numerical and budget refusals leave all accepted owners unchanged.
Transport payloads remain attached to the accepted candidate.

This is a fixed-geometry, held-normal-load phase with constant sampled coefficients
inside its native curve. Endpoint thermal storage and phase-average coefficient
sampling require time refinement, including the timing of heat released before
arrest. Source regularity and physical branch selection remain necessary;
nonlinear convergence does not establish unique source roots or spatial thermal
accuracy. G1/G3/G4 tests cover independent finite-cell balances through arrest and
rest, temperature-dependent stopping, full-tensor/frame-covariant curved travel,
original query receipts and atomic refusal/retry.

## Spatial shock loads on original structural pressure surfaces

The opt-in `experimental-shock-loading` feature supplies `shock_loading`.
Frozen and species-bearing one-dimensional Euler adapters consume the accepted
SSP RK wall momentum exchange, retaining its actual accepted duration. They
require reflecting, stationary, adiabatic ends with zero mass/energy transport.
An explicit global axis fixes orientation; held exterior absolute pressure is
subtracted once from the structural load, never from the gas flux. The caller
supplies the most recent report; a mismatch with current extensive state refuses.
The optional `experimental-compressible` feature adds the original 2-D and 3-D
Cartesian Euler owners' latest per-face wall impulses, with orthonormal frame
embedding. `WallPressureImpulse::from_euler3d` retains all six boundary sides
and both tangential face coordinates. A nonuniform volumetric pressure footprint
therefore retains both transverse moments instead of becoming an averaged load.

Cold surface preparation reuses `fs-plate::ShellPressurePort` and retains sparse
original-node P1 pressure rows. Each explicitly selected planar patch belongs
to one fluid face. Duplicate facets, area/normal mismatch, differing accepted
intervals, work-budget failure and cancellation refuse without a partial load.
Unloaded and supported nodes keep their original identities. Returned impulse,
interval-mean force, resultant, torque and spatial virtual work preserve this
discrete pressure map. Per-face Cartesian loads retain nonuniform shock moments.

This is a stationary-wall load transfer and a one-way prescribed structural
forcing path. It does not move the gas boundary, assert joint fluid/solid energy
closure, resolve within-step load variation, model fracture or validate material
shock response. G1 source-Hugoniot reflection and independent fluid/support
momentum, G0/G3 affine surface work/moment and native modal impulse-consumer
tests accompany explicit boundary/geometry/refusal and cancellation checks.
The 3-D tests independently reconstruct all-six-wall fluid/support momentum,
rotated per-face resultants and both moments, and affine P1 virtual work with
exterior pressure. They also check retained loads after cancelled source steps
and refusal of periodic/open faces, bad embeddings and invalid face indices.

`WallTractionImpulse::from_molecular_euler2d` and `from_molecular_euler3d`
consume the source-owned molecular Euler operators' complete accepted loads,
including both molecular halves and the middle advective step. They never
substitute the borrowed Euler owner's partial report. The gas-side normal
traction is signed `p-tau_nn`; tensile viscous loads remain negative rather
than being clipped or represented as absolute pressure. `WallPressureImpulse`
retains its positive Euler-pressure meaning. Both types share one private
normal-load account and the same original P1 projection kernel, exposed as
`project` for pressure and `project_tractions` for signed traction. No temporary
pressure-vector allocation or extra force distribution is introduced.

The molecular adapters still require stationary adiabatic slip walls with zero
mass, energy and tangential momentum transfer. Explicit exterior pressure is
subtracted over the same accepted interval. Actual transport-law tests exercise
a tensile 2-D wall load through nodal force and affine virtual work, and a full
3-D six-wall momentum balance with retained loads after late cancellation.
Both cases differ from their Euler-only middle-stage impulses. They transfer
normal traction only; moving walls and reciprocal fluid/solid energy closure
remain separate capabilities.

`WallVectorTractionImpulse::{from_molecular_euler2d, from_molecular_euler3d}`
and `ShockLoadSurface::project_vector_tractions` additionally preserve every
component of the complete source-owned wall momentum. The 3-D source's resolved
stationary adiabatic no-slip shear therefore reaches structural drag, torque
and spatial virtual work through scalar P1 area weights prepared once from the
same original surface rows. Normal stress and transverse impulse remain in
their original Cartesian decomposition until exterior pressure is subtracted,
so a large balanced ambient pressure cannot erase small rotated shear. No
traction is inferred from scalar pressure. Compensated shared-node assembly
and resultant, moment and work reduction retain small transmitted impulses
among opposing face loads. Both planar and volumetric sources preserve native
no-slip shear; the planar face area includes its physical depth exactly once.
Authored moving or isothermal walls refuse even at zero net
energy transfer; changing a molecular wall configuration invalidates its prior
complete-step report, while an identical configuration preserves it. Independent
resultant, pure-couple, affine-work, frame-rotation and native no-slip half-cell
shear checks cover this stationary one-way load transfer. It does not advance
structural motion or return that motion to the fluid.

## Constitutive feedback through public material equations

`constitutive_graph::MaterialGraphResidual` closes caller-authored, dimensioned
balances of original law inputs and outputs on the existing Newton--Krylov
solver. Scaled unknowns and residuals keep physical units explicit. Jacobian
actions use the material graph's exact chain rule at the same committed
history as each residual; graph-internal edges remain acyclic. A finite work
budget, latched original source errors, cancellation and independent final
equation tolerances precede returning a complete candidate. No physical
history mutates during a trial or failed solve. The native registered
hyperelastic material-card consumer and a closed two-variable electrothermal
reference exercise this API. Convergence of arbitrary caller balances does
not confer conservation, passivity, calibration or a time-discretization law.

## Additional physical links in the shared material heat assembly

The material, electrical-contact and complete-heat prepared owners expose
`prepare_with_link_conductances` / `resolve_with_link_conductances` paths for
another physical law to supply native thermal links. Ordinary Fourier bindings
and the supplied conductance map must partition all active links exactly at
every trial. Additional links require zero authored conductance and finite
nonnegative resolved G; overlap, missing links and inactive IDs refuse.
Default methods retain complete Fourier coverage and cannot silently omit an
additional owner. Preparation retains source pins without caching extra G.

The same original contact work, cell caloric state and thermal residual apply.
Additional link heat is reciprocal internal exchange, never duplicate source
power. The caller retains the source payload and resolves its law at each
proposed state. `with_integrated_conduction_and_link_conductances` preserves
the existing layer-law correspondence checks; integrated layers can replace
only original Fourier links. A material thermal contact cannot also own an
additional link.

The `fs-airflow` finite-fluid convection consumer composes this path with
prescribed reservoirs and the original circuit/enthalpy solve. Its material
viscosity/Cp substitutions change both finite temperatures and native current
and charge while closing the same energy account. These interfaces add no
mass, momentum, pressure-work or time-integration law; those remain with the
participating physical owners.

## Source-integrated conductivity in composed heat transfer

`PreparedHeatNetwork::with_integrated_conduction` routes the existing material
owner's integrated k(T) layer fluxes through the same circuit, thermal and
observer callbacks. Each selected link retains its original conductivity
property/claim, material card, caloric object, fixed conditions and distances.
The native conductance map and reported link secants both contain the applied
flux. Internal heat still transfers once between the original cells; no extra
source or time integrator is introduced. Immutable correspondence is admitted
once, while the original layer owner checks actual phase, support and geometry
on every trial. A separate thermal-contact owner on the same link refuses,
and topology changes require explicit rebinding. The admitted layer model is
quasi-steady, homogeneous, isotropic and pure-phase; it adds no within-layer
storage, tensor transport or moving phase boundary. Focused synthetic tests
compare the nonlinear flux with an independent linear-k integral and verify
that data-only conductivity changes reach finite heat and electrical history.

The same composed owner admits source-backed cylindrical and spherical layer
links through their retained `LayerGeometry`. Bulk-cell distances remain the
physical trace-to-contact radial widths; equivalent resistance lengths stay
inside the material flux owner. The actual interface area, original source
context and caloric/cell correspondence remain authoritative. Geometry is
prepared once, and every subsequent thermal/electrical candidate uses the
original nonlinear radial flux. Focused `material_heat_network` regressions
cover both geometries and reversed endpoint order, joint finite-store energy,
and refusal of resistance metrics mislabeled as physical lengths or incorrect
side widths. This adds no radial storage, symmetry evolution or source law.

M4.1-W1 extends `unilateral_contact` with explicit physical span measures [m],
per-length stiffness, static/kinetic coefficients and a declared friction ramp
speed (mandatory for nonzero kinetic friction). No 1/n normalization or hidden
1e-3 m/s scale remains in the adapter; constitutive refusals retain owner/field
in a typed DContactError. Loss-port Results propagate. `binding::bound_point_obstacle`
lowers SourcePointContact verbatim and is also used by MaterialContact.
`BoundNormalStorage` composes the existing point or finite-gap potential with
the pHS Storage contract, preserving original U/F/dF and source errors in its
fallible API; the infallible Storage view signals invalid trials as nonfinite
values for the numerical owner's refusal. Whole finite-gap footprints already
integrate their areas and are not multiplied by another area weight. Sources
must have an unloaded zero-approach origin. Point Hunt--Crossley loss is a
separate port evaluated at actual discrete effort and explicit dt; its envelope
is checked before projection. No time integration, cancellation framework,
material data, thermal evolution or interpolant authority is invented.

Multiphysics composition through a versioned port-thermodynamic vocabulary: a
lossless interface relation, explicit storage/dissipation/source primitives,
a complete stream-flux bundle, and evidence-bound closed-window physical
accounting without false passivity.

## Reacting atmosphere on the native moving-liquid timeline

`GasReaction::liquid_candidate` in `liquid_gas` composes the original closed gas
reaction owner with `fs-lbm`'s finite-atmosphere/VOF enthalpy owner. The native
liquid timestep supplies one interval to chemistry, thermal exchange and flow.
The reaction changes composition at conserved total internal energy, retaining
species formation references. Its accepted original gas sources then supply
film conductivity, caloric inversion and pressure during heat exchange. The
same liquid source laws feed native viscosity, buoyancy, conductivity and any
attached current/heat response. Subsequent calls use the actual cooled gas and
remaining species inventory.

The candidate clones native populations, carry, source/contact histories and
finite walls; no caller-owned state changes on source, thermal, reaction,
combined-energy or cancellation refusal. An independent final account measures
gas internal energy, post-advection liquid enthalpy and finite wall enthalpy
against authored/electrical/physical viscous heat. Interfacial heat cancels once;
formation energy is not credited again and numerical kinetic losses are not
thermal sources. Cancellation is checked at existing source and solver stages
and immediately before candidate publication.

This is a first-order reaction/thermal/flow split, with well-mixed fixed-volume
gas, the original constant-density liquid approximation and explicitly exposed
film geometry. It adds no gas traction, flame or shock propagation, evaporation,
phase conversion, moving wall or global mechanical/entropy theorem. Time and
space refinement remain necessary. `tests/reactive_liquid_gas.rs` supplies G1
independent isomer/finite-capacity balances, G3 data-only changes that reach
native VOF populations, continuation and G4 rollback checks; synthetic fixtures
do not constitute empirical combustion or liquid validation.

`GasReaction::liquid_surface_candidate` shares that chemical staging and
whole-interval audit with the native solved liquid/gas/radiation interface.
Product composition, NASA-9 energy and transport feed its surface-temperature
solve; the original liquid `sigma(T_surface)` then drives native capillary
pressure and Marangoni traction while bulk transport retains its enthalpy
temperature. Radiation enters the external account once, and gas/liquid
transfer cancels within the combined finite inventory. The complete-surface
geometry and source restrictions, fixed gas mass/volume and first-order split
remain binding; this adds no surface-energy or mechanical-energy theorem.
`tests/reactive_liquid_gas/surface.rs` supplies independent zero-storage
radiative and finite-capacity references, data-only kinetics and surface-tension
changes reaching native populations, continuation and refusal rollback.
These are manufactured model checks, not experimental validation.

`GasReaction::pressurized_liquid_candidate` and
`pressurized_liquid_surface_candidate` additionally select the native finite-gas
pressure owner. The EOS pressure after chemistry and thermal exchange drives
the original gas-link reconstruction during exactly one liquid flow step.
They require an explicit fixed pressure reference, a bounded lattice-density
excursion, stationary solid geometry, and one well-mixed atmosphere connected
to every native gas link. Gas volume is a physical input; numerical fill never
supplies compression work. The existing thermal audit remains distinct from a
complete mechanical-energy or entropy theorem. G1/G3/G4 tests in
`tests/reactive_liquid_gas/pressure.rs` check independent chemical/film pressure
and native momentum, radiative feedback, the single clock and refusal rollback.

All reacting thermal entries preserve the native wall source accounts. Gas may
exchange heat only through exposed finite-solid faces; its subsequent contact
heat reaches liquid without a fabricated direct gas/liquid film. Radiation
deposited in those solids contributes once to the enclosing external supply.
Gas/solid and solid/liquid transfers remain internal. The finite-wall regression
checks both direct and wall-only gas paths with sourced wall radiation.

The enclosing gas storage change is the sum of the original chemistry energy
increment and the original finite-atmosphere thermal increment. It is not
recomputed by subtracting absolute formation-referenced endpoint energies.
The frozen-composition thermal path therefore retains the native owner's
reference-invariant heat account, including large common energy-datum shifts.
This does not extend the changed-composition reaction owner's numerical scope.

### Finite-rate condensation and evaporation on native moving liquid

`liquid_gas::evaporation::advance_liquid_with_gas` joins the original
`fs-material::evaporation::gas` Stefan-film/caloric owner to actual native VOF
mass and enthalpy. One caller-measured Cartesian patch binds an occupied liquid
interface to adjacent in-domain gas. The original native caloric chart supplies
temperature; its tracked mass selects a reproducible transfer before the gas
caloric solve. The SAME selected mass debits the finite gas and credits native
liquid. Latent and upstream-vapor sensible heat enter once, using the actual
native post-transfer mass. No second liquid-film inventory is retained.

The original native source refresh and momentum/enthalpy transport then advance
one retained timestep. The immutable candidate publishes both participants
together. Actual total and volatile mass balances include subsequent native
transport drift; the source-relative thermal defect includes phase caloric
reconstruction, native caloric representation and native advection separately.
Original gas-species rounding remains separately mass-budgeted. Native external
heat stays in its existing flow report and is not mistaken for latent heat.
Attached electrical cores require the original conduction timeline and refuse
this direct flow entry point. `tests/evaporating_liquid_gas.rs` covers condensation warming,
evaporation cooling, finite backpressure on continuation, independent whole
inventory energy/mass checks, reference shifts, capacity/geometry refusal and
atomic cancellation. Its authored phase and transport data are manufactured
coupling checks, not empirical water-rate validation.

`advance_liquid_with_gas_and_heat` additionally uses the original implicit
reciprocal gas-film heat solve before phase transfer over the SAME interval.
Its explicit sensible coefficient excludes vapor-carried sensible and latent
heat. The caller supplies the original pure-liquid enthalpy branch; temperature
still comes from the native caloric chart. Native mass selection occurs only
for the actual phase debit after heating. The final native enthalpy includes
both heat and phase, followed by exactly one native flow step. Reports retain
the separate heat receipt, phase-only receipt after heat, and original-to-final
`exchange` receipt, with one combined energy budget. A zero heat coefficient
retains the original phase-only endpoint and bypasses unused heat controls.

This bounded first-order split needs timestep refinement. Gas occupies a fixed
connected control volume; native liquid volume follows its original constant
density. Phase-law recession is diagnostic, not a second interface update.
The phase-only entry point supplies no separate sensible film conduction.
Neither path infers sealed-vessel displacement, gas traction/recoil,
nucleation or phase-front overflow/exhaustion remap. Existing
solved-surface stress and dynamic-contact histories retain their native refusal
boundary. Exposure is held during this one interval and re-admitted on the next.

`advance_liquid_with_pressurized_gas` and its `_and_heat` variant, plus the
corresponding `advance_liquid_surfaces_with_pressurized_gas` routes, pass the
actual post-exchange finite-gas EOS pressure to the original native
`step_with_heat_and_gas_pressure` owner during that same flow step. The original
gas pressure receipt records the phase-initial pressure, applied endpoint
pressure, fixed physical volume and explicit density-excursion allowance.
All exposed links must belong to the same caller-declared well-mixed atmosphere.
A fixed pressure reference, finite allowance in (0,1), stationary solid geometry
and the original per-link density/pressure gates remain required. Pressure is
mechanical loading, never an extra heat debit; sealed-vessel displacement,
compression work, recoil and a global mechanical-energy theorem remain outside
this fixed-volume split. Late pressure refusal discards both phase endpoints.

`liquid_gas::evaporation::reactive` connects `GasReaction` to the single- and
multiple-surface phase owners, each with an optional reciprocal sensible-heat
route. Chemistry uses the same retained native timestep at conserved total
internal energy. Its actual composition and temperature feed phase exchange
from the identical original source catalog; native flow advances once.
Formation energy is never a second heat source, and chemical production of the
volatile species is distinct from its phase transfer. Original source, donor,
geometry, mass-representability and attached-owner limits remain active.

Chemical-intermediate and final temperature bounds plus sums of absolute
chemical/phase/native mass and energy defects gate atomic publication.
Cancellation or refusal preserves every accepted owner. This is a first-order
reaction/heat/phase/flow split with a well-mixed fixed-volume gas. The original
unpressurized entries add no gas traction, sealed-vessel displacement, flame or
global entropy theorem. Five
G1/G3/G4 regressions in `tests/evaporating_liquid_gas/reactive.rs` cover chemical
reversal of phase direction, independent whole inventories, all four routes,
original-owner agreement, continuing species, equilibrium/zero-heat limits,
canonical surface ordering and rollback.

The four `GasReaction::pressurized_evaporating_liquid*_candidate` entries compose
the same chemical transaction with the original pressurized single/multiple
phase owners, including their sensible-heat variants. Chemistry changes the
actual composition before phase transfer; pressure uses the final mass and
temperature after every exchange. The phase pressure receipt starts after
chemistry, and one native interval advances. Native pressure and reactive
pressure regressions check original-owner equivalence, a pressure-response
oracle, opposing surface transfers, canonical order, conservation and exact
retry after late pressure rejection. These routes retain the same fixed-volume
and first-order physical limits; numerical regressions require native execution
and do not establish empirical constitutive validation.

### Native phase transfer through full finite-gas conduction

`liquid_gas::evaporation::conduction::advance_liquid_surfaces_with_gas_and_conduction`
connects the original shared-gas phase transaction to the complete native
thermal/flow owner. The represented phase mass first updates the same native
liquid caloric field; the exact post-phase closed parcel then enters its
original isochoric finite-gas owner without reconstruction. Bulk liquid
conduction, finite solid storage and contacts, gas films, wall radiation,
attached electrical current, and optional physical viscous heat all remain
on that original timeline. There is no duplicate sensible-film coefficient
or second native flow step.

The candidate retains liquid, gas, finite solids, original catalog and both
original reports. The phase report already contains the native thermal/flow
defect; the complete energy account adds only the original gas caloric defect.
Absolute local defects share the whole interval budget. The thermal stage
preserves gas mass, composition sources and physical volume. Optional pressure
uses the actual post-thermal EOS; its receipt starts after phase transfer.
Every late source, work, thermal, pressure, budget or cancellation refusal
discards all candidate owners together.

`GasReaction::conducting_evaporating_liquid_candidate` adds the existing
chemical transaction to that same connection. Chemical formation energy stays
with its original owner, chemical production of vapor remains distinct from
phase transfer, and native flow advances once over the retained interval.
The enclosing chemical account consumes the lower complete mass/energy
defects. Chemical, post-phase and final gas temperatures each obey the declared
departure bound. Finite walls publish in the same reactive candidate.

These are first-order phase/conduction/flow and
chemistry/phase/conduction/flow splits requiring complete-timestep refinement.
The original fixed-volume atmosphere, stationary geometry, source, donor,
caloric, pressure-reference and attached-history restrictions remain binding.
This connection supplies no compression work, sealed-vessel displacement,
dry-out/remap closure or complete mechanical-energy/entropy theorem.
Regressions in `tests/evaporating_liquid_gas/conduction.rs` and
`reactive/conduction.rs` cover opposing phase streams, finite wall/gas/radiant
accounts, a real attached electrical core, post-thermal pressure, original-owner
agreement, independent whole inventories, one clock and late refusal/retry.
They require native execution before any numerical validation claim.

### Retained moving gas through native thermal and chemical owners

`liquid_gas::evaporation::momentum::conduction` connects distributed phase
transport to the same finite-gas/wall conduction, liquid and wall radiation,
attached electrical heat, optional viscous heat and native pressure/flow owner.
Every phase donor uses its actual retained velocity. Receiver mixing heats the
original receiver before thermal exchange; the original post-phase gas enters
the thermal owner without a new EOS or caloric model. Exact gas mass, volume and
composition checks at the thermal boundary make preservation of its post-phase
bulk momentum also preserve its kinetic inventory. Optional pressure uses the
actual post-thermal EOS in the one native flow step.

`GasReaction::moving_gas_evaporating_liquid_candidate`, the distributed and
sensible-heat variants, and the conducting/radiative-conducting variants use
the original reaction integrator and shared reactive admission path. Homogeneous
chemistry preserves exact original parcel mass and applies no external impulse.
The moving adapter explicitly verifies mass and volume before reattaching the
same momentum; formation energy remains in the original chemical/internal-energy
owner. The reacted donor inventory then enters the original moving phase owner.
Chemical, post-phase/mixing and final thermal gas endpoints are checked against
the original reactive temperature allowance. All finite walls publish with the
liquid, gas and its momentum only after the complete transaction succeeds.

Signed and absolute defects retain original accounting boundaries: the single
phase receipt contributes its mechanical defect once, distributed phase already
includes it, and complete conduction adds only the additional finite-gas caloric
defect. There is one physical interval and one native flow step, with no second
chemical heat, receiver-mixing heat or whole-thermal-flow defect. Existing
canonical surface order, deterministic original solvers, source/geometry/work
admission and cancellation checkpoints remain binding.

These first-order fixed-volume compositions retain a source-only momentum
account. They do not infer resolved gas flow, wall drag, photon momentum,
compression work or a closed mechanical balance for native pressure/wall forces.
G1/G3/G4 fixture definitions in `conduction/momentum.rs` and
`reactive/momentum.rs` specify opposing phase streams, independent caloric plus
source-kinetic inventories, exact original-owner chain agreement, one clock,
post-thermal pressure and atomic late refusal/cancellation. Those full material
integration fixtures were authored but not compiled or executed in the delivery
environment; parsing/source review is not numerical validation.

### Liquid radiation on complete reactive and phase timelines

`advance_liquid_surfaces_with_gas_and_radiative_conduction` passes the original
resolved-cell-temperature liquid radiation law to the same finite-atmosphere
thermal owner. `GasReaction::radiative_liquid_candidate` and
`radiative_conducting_evaporating_liquid_candidate` retain that boundary through
the original chemical and phase accounts. The original nonradiative entries
share the implementation with an empty radiation slice.

Liquid and wall radiation are external supply exactly once; finite gas films
and solid/liquid contacts stay internal. A radiation-only liquid boundary can
leave the gas's thermal storage unchanged. Original emissivity selection,
surface geometry, combined thermal stability and optional endpoint pressure
remain authoritative. The caller must resolve the liquid thermal boundary
layer or justify cell temperature as its radiating-surface temperature.
This does not supply the distinct interface-temperature/latent-power equation
needed for solved-interface phase exchange; attached solved-surface and wetting
histories retain their existing refusal.

### Original gas plumbing on the native phase/thermal clock

`evaporation::transport::advance_liquid_surfaces_with_gas_transport` first
advances the original finite-gas transport snapshot over the liquid's retained
timestep. It passes one explicitly nominated atmosphere node's exact transported
parcel into the original phase and radiative-conduction transaction. The complete
gas map, native liquid, finite walls, source catalog and both original reports
publish together. Idle gas nodes remain in the result.

Every original transport path shares its original gross-withdrawal allowance.
The final gas map must pass the original endpoint driving-force and source
checks after phase and heat exchange; an active pressure or permeation reversal
rejects the entire candidate. Original source-relative transport residuals and
absolute local defects join the complete phase/thermal budgets. Transported
enthalpy, phase latent heat and internal film/contact heat are not reapplied
as external sources.

This is first-order transport/phase/conduction/flow splitting. Arrivals can
participate in the later phase stage but cannot fund another outgoing path
within the same original transport snapshot. It adds no simultaneous
phase/plumbing solve, jet momentum, spatial gas CFD, compression work or
whole-system entropy theorem. Regression sources cover real-nozzle phase
feedback, independent network/liquid/wall inventories, actual final-pressure
refusal, empty-path equivalence, retained idle nodes and deterministic retry.
Native execution remains necessary before numerical validation claims.

## Assembled shell adhesive dissipation

`shear_shell::assembly::MaterialShellAssembly::reduce` projects the original
mortar glue damping through the same assembled mass-normalized modes as
stiffness, component viscosity, bridge forces and pressure surfaces. It sums
the complete modal tensor before the existing positive-semidefinite damping
owner realizes passive channels, preserving cross-mode terms and undamped
subspaces. With fixed stiffness and geometry, changing glue viscosity changes
neither the elastic basis nor the radiation map. Zero supplied viscosity adds
no channel; no loss is inferred from a material name.

`has_viscous_glue` reports supplied positive edge or face viscosity. Projection
and connection budgets include this work. Only rigid mortar ties eliminate
slave coordinates; supports on compliant slaves remain active. Original tie
owners still validate coefficients and geometry. Source adapters own the
validity of a constant-coefficient glue law across the retained state and
frequency range.

## Purpose and layer

Layer L3 (multiphysics coupling). Domain assembly depends on physical
solvers. The neutral vocabulary is owned by `fs-ports` and reexported here
without conversion or duplicate types. That lower boundary depends only on
`fs-qty` and `fs-iface`; FEEC imports it without depending on domain assembly.

## Public types and semantics

### Shared cavity and distributed downstream duct

`SharedCavity::with_downstream_duct` binds the unchanged compact cavity/neck
to an additional source-bound `DgPortLoad`. Admission requires an
`ExteriorMouth` neck with no direct exterior load, equal inlet gas and area
(roundoff-only area allowance), an ordinary non-driven inlet, cold accepted
clocks and identical step duration. The original neck retains its bore, inner
end inertia and selected viscous law. DG owns additional distributed length,
wall memory, branches and outlet terminals.

A scalar displacement mobility is prepared once. Each step solves mean inlet
pressure against current complete DG history; the original modal and DG
candidates publish together only after finite pressure, swept-volume,
equal/opposite port work and joint budget checks. The original energy ceiling
covers both owners' represented storage; the pressure ceiling also bounds mean
DG inlet pressure. Rational terminal work remains signed boundary exchange.
Dense-DG loss inferred from its storage deficit is explicitly reported as such,
not as an independently integrated heat channel.

The held-pressure/held-mean-flow connection approximates the continuous
interconnection; it is not the joint continuous exponential. Coupling-clock,
spatial and constitutive-band convergence remain separate. Tests replay both
original propagators with accepted interface quantities, retain actual viscous
and thermal memories, refine the coupling clock at fixed geometry, and refuse
late candidates without publishing either owner. The returned cavity accessor's
harmonic response describes the cavity subsystem and excludes the attached DG
load. A broadband 5 cm/order-2 fixture had nonmonotone 24-to-48 kHz
endpoint error. Its retained DG spectrum reached 47,649 rad/s; the unchanged
error gates passed at 96/192/384 kHz against a 1.536 MHz fixed-spatial reference.
This is fixture convergence evidence, not a universal timestep certificate or
a continuum validation.

### Molecular exterior loading of original structural motion

`BemLoadedModes::prepare_molecular_harmonic` binds the model's original solid
or finite-thickness shell skin to `MolecularExterior` and optional original
`LocalBoundary` laws. The equilibrium gas density and sound speed must match
the source model. Every original mass-normalized structural coordinate, full
modal viscosity and Maxwell memory branch is retained. Driven panels remain
explicitly rigid under the wall owner's admission; other patches may absorb.

The complete unit-source boundary bank gives the discrete acoustic load Z and
matched finite-distance pressure/particle-velocity transfers. The mechanical
block is `K - omega^2 I - i*omega*(C + Z)` under `exp(-i*omega*t)`, with the
unchanged Maxwell rows appended. Only the boundary variables are eliminated;
this is the exact Schur complement of the selected discrete source equations,
not a fitted impedance or additional source-mode truncation. The original
boundary bank and one balanced mechanical/memory LU are reused across all
admitted coherent force cases. No independent receiver solve or absorption
filter is added.

`PreparedMolecularLoadedModes::harmonic_response_batch` returns force-scaled
displacement, velocity, complete elastic memory, fluid reaction `-Z*v`, matched
receiver fields, and original unscaled mechanical/memory residuals. Cycle work
separately reports input, full Kelvin heat, Maxwell heat, exterior work and
wall absorption. Exterior work includes eventual molecular heat and outgoing
flux; it is not distant acoustic flux alone. All cross terms and raw signed
channels remain, with no passivity repair or duplicate gas-loss debit.

The source owns geometry and mechanics throughout preparation and replay.
Positive frequency, the original material band/formulation, nonempty physical
receivers, and separate boundary, source-bank, mechanical scratch, output and
work budgets remain explicit. Cancellation or a refused force case publishes
no partial response and leaves the cached preparation reusable. This is a
stationary homogeneous molecular harmonic response, not a time-domain material
realization, complex resonance continuation, mean heating or boundary-layer
model. Original spatial discretization limits still apply.

Five focused regressions cover an actual FE-derived shell with full damping
and material memory, loaded-motion and receiver agreement with the existing
full pressure-state solve, coherent forcing and passive wall work, a physical
mass/Kelvin/Maxwell breathing sphere on three meshes with an independent
complex-wave analytic solution, empty-reservoir recovery, and budget/source/
cancellation refusal. These test definitions are not native execution evidence;
native validation remains pending the shared dependency/resource blocker.

`PreparedMolecularLoadedModes::harmonic_incident_response_batch` extends that
same source to incoming plane waves, finite monopoles, force dipoles and mixed
coherent fields. It projects the stationary illuminated body's total pressure,
including diffraction, into blocked acoustic force on every original mode.
This force and the supplied external mechanical force share the cached loaded
mechanical/memory solve. One additional blocked-body BEM assembly/LU serves the
incident batch; the earlier motion boundary bank is reused without another
motion-field solve. The original gas, wall laws, formulation, source modes and
receiver declarations remain fixed by preparation.

The response exposes actual motion, independently projected total fluid
traction, total boundary pressure/normal velocity, and total receiver pressure/
particle velocity. Blocked, motion, interference and total exterior/wall work
are separately contracted from their physical traces. In particular, total
outward boundary work is signed under illumination: sound-only excitation can
make it negative while structural or wall heat is positive. The balance is
external mechanical input = Kelvin heat + Maxwell heat + total exterior work
+ total wall absorption. Incident/motion cross work is neither omitted nor
recovered by subtracting large powers; the existing motion-only response and
its effective-force ledger remain explicitly distinct from these total fields.

`MolecularIncidentLimits` adds source clearance, analytic source-point work,
projection/observation work and conservative live-payload caps. These and the
inherited force-batch limits run before the blocked solve. Failed batches
publish no partial response and leave the prepared owner reusable. G0/G1/G3
regressions cover finite-source coherence with original Robin walls, direct
combined-drive field parity, exact acoustic/mechanical force cancellation,
resource/cancellation refusal, and an independent three-mesh compliant-sphere
oracle for blocked force, material-memory motion, paired fields and signed
total work. They are authored but remain unexecuted natively under the same
dependency/resource blocker; no new continuum or performance claim is made.

`BemLoadedModes::prepare_molecular_harmonic_adaptive` selects the original
adaptive retained-triangle pressure/velocity receiver for the same loaded
mechanics. A single receiver row bank serves every structural source column;
subsequent force cases reuse the matched transfer bank. Each response retains
absolute quadrature estimates for its actual coherent modal velocity. This
selection changes neither the boundary load nor the original mass, stiffness,
Kelvin/Maxwell state or surface-work calculation.

Sound-driven responses retain that same receiver choice. One additional
source-bound receiver bank serves the entire blocked incident batch, while
the existing motion receiver bank remains cached. Total pressure/velocity add
blocked and motion fields coherently; total quadrature estimates add the
blocked estimates to absolute-rate-weighted motion estimates. These estimates
do not cover BEM discretization, material/source error or relative accuracy
near a field null. The new row-storage and quadrature-work limits remain
separate. Adaptive force outputs reserve another 32 bytes per case/receiver;
adaptive incident payload replaces `cases*8*r+4*r` with `cases*12*r+6*r` in the
existing complex-entry reservation. Original centroid reservations are unchanged.

A focused G1/G3 regression uses original two-mode mechanics, Maxwell memory
and local wall physics with multiple force/illumination cases. It checks
unchanged loading, motion and work, direct original-owner adaptive field
agreement, shared-row reuse and total discrepancy addition. The test is
authored and syntax checked; native execution remains pending the existing
dependency/resource blocker.

### Open-shell forced response with two-sided diffraction reaction

`bem_loaded_modes::open_screen::OpenScreenLoadedModes::from_shell` binds the
original `ShellReduction` midsurface and every mass-normalized translation to
the native open P1 pressure-jump owner. It retains the complete zero-state
potential Hessian K and the supplied admitted full `ModalPlateViscosity` C;
matching viscosity dimensions does not itself establish source-basis provenance.
Both faces contact the same homogeneous exterior. The selected acoustic chart
has zero thickness and no cavity-backed face. Physical source thickness remains
in the mechanical mass/K/C, and rotations have zero midsurface offset.

`prepare` preserves the original open topology, dense-screen ceiling and
independent edge/triangle/angular resolution and work admission. Optional
conforming acoustic midpoint subdivisions interpolate the source's facet-local
P1 motion exactly. Adjacent folded normals are not averaged and free-edge drive
velocities are not zeroed. This refinement leaves the source's piecewise-flat
geometry and structural basis fixed; it supplies neither edge grading nor a
structural convergence receipt.

`PreparedOpenScreenModes::harmonic_response_batch` solves the real-frequency
equations `D = K - omega^2 I - i*omega*(C + Z_screen)` in `exp(-i*omega*t)`.
The original boundary assembly/LU serves all source columns, and one original
cancellable complex LU of D serves every complete authored force case. All
self/mutual entries and source signs remain in place. Fluid reaction is
`-Z_screen*modal_velocity`; no area conversion, duplicate wet face, fitted
damping or hidden broadening is applied. Returned force-scaled displacement,
velocity, reaction, source structural energy and original mechanical-row
residual accompany a separate numerical balance of external input, complete
viscous heat and exact P1 pressure-jump work. Raw signed boundary work is not
clipped, symmetrized or presented as a continuum passivity certificate.

`OpenScreenLoadedModes::from_shell_molecular` derives the equilibrium medium
from the original evaluated gas. The prepared owner's
`harmonic_molecular_response_batch` binds that same equilibrium density/speed
to the complete molecular compression law at the physical frequency. It uses
the original source K/C, exact P1 work map and shared mechanical solver; only
the source-owned boundary matrix changes to H(k(omega)). Every source column
shares one boundary LU and every complete force case shares one mechanical LU.
Self/mutual terms, free-edge velocity and the original modal basis remain intact.
`MolecularOpenScreenHarmonicResponse` reports actual motion, reaction, structural
energy and original row residuals, with the signed pressure-jump input in
`power.exterior_w`. That includes molecular gas heat and outgoing energy once;
separate sheet loss and structural Maxwell memory are absent from this screen
owner and their power fields are zero. No far-field flux is inferred from it.

`MolecularOpenScreenHarmonicBatch::prepare_receivers` retains one adaptive
complex P1 row bank borrowed from that original solved source. Every case uses
its actual loaded modal velocity for coherent pressure and all three particle
velocity components on either side, without another boundary solve or receiver
quadrature. Independent receiver storage/work/depth/clearance limits remain
explicit. Absolute estimates describe integration of the fixed solved jumps,
not mesh/material uncertainty or cancellation-relative accuracy. The original
full mechanical batch limits still apply before boundary work begins. Two
focused regressions compare native empty-reservoir response, independent 2x2
mechanical inversion, complete mutual reaction, coherent fields, signed power,
source mismatch and cancellation. They are authored and syntax checked; native
execution remains pending the existing missing dependency/build-runner blocker.
This frozen harmonic model does not add cavity backing, finite acoustic
thickness, shear/thermal boundary layers, a time-domain realization or molecular
loaded-pole discovery.

`PreparedOpenScreenModes::harmonic_molecular_incident_response_batch` drives
the original shell with incoming sound and an optional authored mechanical
force (an explicit zero vector for purely sound-driven cases). It uses the
original blocked pressure jump's work-conjugate force, then solves the unchanged
complete K/C and molecular self/mutual reaction. One boundary H/LU serves all
motion columns and blocked illuminations, and one mechanical LU serves all
resulting force cases. Original source phases, geometry and material remain
bound to the returned batch; no external impedance matrix is substituted.

`MolecularOpenScreenIncidentPower` independently contracts external mechanical
input, blocked acoustic input, original structural viscosity, motion-only
exterior work and TOTAL signed pressure-jump work. The original total fluid
reaction is blocked force minus Z times actual velocity. Thus external input
equals structural heat plus total exterior work; a sound-only damped shell has
negative total exterior work. The blocked pressure/motion exchange is directly
contracted, never recovered as a difference of large radiation powers. This
zero-thickness impermeable screen has no separate wall-absorption channel.
Normalization, compensated work sums and explicit physical balance refusal
preserve small channels without clipping signs or claiming continuum passivity.

Prepared incident receivers use the actual sound-driven modal velocity with
the original incident + blocked-scattered + motion field. One shared P1 row
bank serves every case; the analytic incident field is cached and added once.
The full mechanical response, force-column/conversion payload and additional
work are admitted before boundary assembly. Two focused regressions exercise
sound-only and joint force/sound motion, signed work, independent existing-owner
mechanical recovery, coherent mixed-source receiver fields and transactional
resource/cancellation refusals. They are authored and syntax checked; native
execution remains pending the existing dependency/build-runner blocker.

`OpenScreenHarmonicBatch` retains the original unit-source jump solutions.
`pressure_jump`, `pressure_at`, `far_field` and `radiated_power_w` combine them
with the actual loaded modal velocities before the original field integration,
preserving coherent phase and cross terms without another solve. Adaptive
near-field receiver/clearance/kernel-work limits and the original angular rule
still apply on both sides. The independently integrated angular power remains
distinct from pressure work; the maximum unit-source backward residual is
neither a coherent cancellation error bound nor a spatial error estimate.

Source modes are capped at 256 and never truncated. Before a frequency solve,
the declared mechanical live-payload allowance includes
`cases*(48*modes+512) + 48*modes^2 + 160*modes + 4096` bytes for responses,
original/scaled matrices, LU and scratch. Boundary/source storage is separately
bounded; these logical counts exclude allocator overhead and process RSS.
Mechanical work is bounded by `32*modes^3 + 256*(cases+1)*modes^2`, with all
native boundary work additional. Source tangent assembly uses the original
cancellable Hessian loop. A failed batch returns no partial result.

Four authored physical-source regressions compare full reaction with an
independent two-coordinate dynamic-stiffness inversion, expose omitted mutual
loading, verify coherent jumps/two-sided fields/power, and pin work/payload,
native topology/band and cancellation/replay refusals. They have not run under
the current missing canonical dependency checkout. This is a linear forced
response about the reference source; it does not provide nonlinear motion,
a time-domain screen realization, a finite-thickness or
interior-cavity model, or the independent spherical/open-cylinder acceptance
required by M2.5-W4b.

`PreparedOpenScreenModes::resonance_pencil` adds loaded-mode analysis of the
same complete K/C source and real-frequency pressure-jump discretization. With
G in mode-major order, its original operator on `[x, jump_unknowns]` is
`[K-omega^2 I-i omega C, G; -rho omega^2 G^T, H(omega/c)]`.
G is the exact existing facet-local P1 pressure/velocity work pairing, prepared
once. The analytic derivative keeps every block, including `H'(omega/c)/c`.
No pressure unknown or zero-opening structural coordinate is removed. A contour
therefore sees the holomorphic augmented operator, not the meromorphic quotient
`H^-1 G^T`. Local initialization alone solves H for a supplied displacement
seed and refuses a singular isolated boundary block.

`OpenScreenResonancePencil::{original_pencil,resonance_pole,resonance_spectrum}`
reuse the existing bordered Newton and Beyn owners. Fixed reference balancing,
physical eigenstate normalization and compensated original-row acceptance are
shared with the closed-body loaded-mode owner. Complex frequency never changes
the balance. Original source dimensions, callback/seed payload, complete cold
projection and per-callback work, total operator evaluations and full-nodal
output all have explicit limits; native eigensolver scratch/work remains
separately bounded. The whole ellipse and its interior are admitted before
iteration against the boundary's original absolute-wavenumber band and centered
exponential-growth chart. A finite-node contour does not bypass that domain.

Returned modes contain original modal displacement, the full nodal pressure
jump with exact free-edge zeros, mechanical/boundary residuals and native
refinement diagnostics. Amplitude is arbitrary and normalized in the original
coordinates. `decay_rate_per_s` is `-Im(omega)` under `exp(-i omega t)`;
`quality_factor` exists only for a decaying oscillatory pole. Growing accepted
roots refuse. Contour results retain the original rank/count, nested-moment,
conditioning, unrefined-value and evaluation reports. These are numerical
finite-operator results, not interval completeness or physical-root certificates.

The open boundary's cosine-plus-Gram continuation is analytic at zero, so no
closed Burton--Miller `1/k` exclusion is invented. A multiple rigid zero pole
may still refuse native simple-root refinement. This bounded continuation
matches the real-frequency operator; it does not replace the distinct damped
point-kernel Laplace owner, certify all-RHP causality/passivity, supply a CQ
reaction runtime or add material-memory states. Source, mesh, triangle and
angular convergence and physical root identification remain required.
Five authored regressions check complete real-axis harmonic-equation parity,
complex-direction derivatives and zero, an actual non-radiating planar shell
mode's exact decay on three native contour grids, the radiating shell's
independent discrete incompressible added-mass limit, and domain/resource/
cancellation refusals. Package execution remains blocked by the missing
canonical dependency checkout; spherical-shell and open-cylinder validation
remain outstanding.

### Causal modal observation of affine tone envelopes

`BroadbandModalObserver::prepare_affine_tones` and `step_affine_tones`
extend the original continuous-tone observation owner with linear force
envelopes. Both coefficient endpoints share the interval-start carrier basis.
The shared scalar mechanical owner supplies actual motion, input work and
positive heat; the original continuous rational filter receives that entire
trajectory through two additional exact source columns for each unique
frequency. It does not hold an endpoint velocity or fade an already-rendered
pressure. Startup, release ringdown, nearby-tone interference and existing
sampled-delay history remain physical state.

The same immutable source identity, clock, Nyquist guard, fit/held-out gates,
matrix-entry cap, original mechanical limits and receiver-pressure ceiling
remain binding. Constant-tone preparation retains its original smaller path;
a mismatched preparation/step family refuses. Endpoint duplicates are summed
before differencing, matching the mechanical source. A real constant
feedthrough is admitted; imaginary constant quadrature still refuses.
Mechanical state and frame, filter state and delay history commit together
only after all previews and physical gates pass. Refusal changes no accepted
history and never advances caller-owned phase. Observation adds no storage,
heat or reaction load, and its numerical transfer evidence retains the
original estimate-only scope.

### Continuous tone forcing of coupled modal mechanics

`coupled::tones::CoupledModalTones` advances an existing complete coupled
network under simultaneous continuous cosine/sine generalized forces. Force
rows retain the original modal order and N/sqrt(kg) units; coefficients refer
to the start of each accepted sample. Its read-only owning wrapper binds the
prepared operator to the original source, clock and histories. Recovering the
network for a coefficient change requires preparing it again.

The existing exact source assembler retains complete signed spring and
dashpot columns, the prepared coupled-basis direction of constant-loss
dashpots, Maxwell memories, and reciprocal acoustic coordinates. The shared
`fs_phs::linear_exact::PreparedExactTones` integrates this full trajectory,
original source work, component heat, every link's heat, material heat and
radiation heat. Each original link reports its actual negative elastic-storage
change minus its own integrated heat. Spring rest loads enter as DC and their
work is removed from the authored-source account. All coherent cross-tone
terms remain; exactly equal sources combine in caller order with compensated
addition, while nearby tones retain their beats.

The existing state, energy, pressure, material-domain, whole-interval force
and port gates precede the original single publication boundary. The force
gate uses a bound on the continuous forcing amplitude. Refusal preserves
accepted mechanics, memories, radiation, diagnostics, clock and ports.
Real instantaneous pressure-per-velocity maps remain usable; nonzero imaginary
natural-mode quadrature maps refuse. Frequency-resolved broadband pressure,
nonlinear contact and general descriptor evolution remain separate owners.

The complete original and source-augmented coordinates must fit the shared
dense capacity and explicit dimension/byte/work allowance. Original source
setup and physical limits remain binding; nothing is dropped to fit them.
Cold numerical preparation is cancellable. Warm propagation reuses its
matrices and evaluates no exponentials; original transactional staging and
port reporting retain their bounded vector allocations. This finite linear
source path does not establish material, radiation-fit or continuum accuracy.


### Affine tone envelopes through complete coupled mechanics

`CoupledModalTones::new_affine` prepares the original coupled source for
piecewise-affine cosine/sine coefficient envelopes. `step_affine` and its
cancellable variant accept complete mode-major start/end arrays in the same
interval-start carrier basis. The existing phase helper rotates endpoint
coefficients for continuation. Constant preparation keeps its smaller path;
using the wrong step family refuses.

The same source assembler retains original springs/rest lengths, dashpots,
Maxwell branches and reciprocal acoustic coordinates. Start and end duplicates
coalesce independently; rest loads are identical at both ends and their work
is removed from authored forcing as before. Each original heat channel and
signed force work uses the continuous envelope trajectory. The force guard
uses the maximum endpoint coefficient norm per distinct tone, a conservative
bound throughout its affine segment. Additional endpoint buffers and the
four-coordinate source generator count against explicit byte/dimension/work
allowances before preparation.

`BroadbandCoupledModalTones::step_affine` observes this source through the
existing exact observer and shared fitted transfer/delay. It uses the same
stage/material check/observer preview/pressure gate/publication boundary.
Release preserves physical and acoustic histories; it does not zero the
output or reset resonant state. Refusal retains accepted mechanics, work,
heat, ports, filter history, delay and last pressure.

Focused native checks compare a two-body attack/release and its signed rational
receiver against an independently written RK4 force/filter/power ODE; verify
complete material/radiation/source work and heat under subdivision; preserve
bitwise physical accounts when adding observation; and check rollback/retry
and explicit affine-family/capacity admission. This extends prescribed force
reconstruction, not a chirp or a constitutive excitation model.

### Causal broadband pressure from complete coupled tone motion

`coupled::tones::broadband::BroadbandCoupledModalTones` owns one causal
receiver over the complete `CoupledModalTones` trajectory. Its solver-neutral
`ModalTransferSamples` use the original flattened mechanical-mode order,
signed velocity normalization and exp(-i*omega*t) convention. Maxwell memory,
signed inter-component couplings and reciprocal acoustic states remain in the
source operator and change those velocities; they are not added as duplicate
radiators. Preparation requires the source's initial sample and exact clock.

The shared modal transfer fitter now uses `fs_vfit::vector_fit_transfer`,
which permits signed feedthrough while retaining the stable-pole and residual
gates. The impedance fitter retains its separate nonnegative D/E policy.
The receiver realizes each proper rational transfer with its actual real
state space. `PreparedToneObserver` propagates these filter states with the
complete continuous source and source phases. It does not hold an endpoint
velocity, use natural-mode quadrature at a different forcing frequency, or
reconstruct mechanics as uncoupled oscillators. A real constant transfer is
admitted; a complex constant quadrature refuses.

The same shared sampled-delay owner and independent held-out transfer gate
used by the original modal receiver remain binding. Integer delay preserves
sample history exactly; fractional delay retains its explicit interpolation
error axis. Filter and delay histories start at zero, so startup includes the
causal filter transient. `transfer_at` reports the fitted velocity-to-pressure
operator including that actual sampled delay, not a continuum guarantee.

Preparation explicitly clears the original component observation maps on
accepted and candidate source models, preserving all physical coefficients,
kernel identity, load reactions and histories. The new receiver owns the
original network pressure ceiling and returns pressure separately from the
mechanical frame. It contributes no storage, heat or reciprocal force. Source
staging, filter preview and delay preview precede one accepted publication;
any mechanical, material, pressure or cancellation refusal preserves accepted
mechanics, ledgers, ports, filter history, delay history, clock and last
pressure together.

An additional explicit byte allowance covers conservative fitting workspace,
retained rational models, real matrix assembly, filter state and delay
history. The complete augmented source/filter operator must also satisfy the
shared dimension and multiplication-work limits. No mode, tone, pole or
physical memory is removed to fit these bounds. Warm filter propagation uses
cached matrices; the original source's transactional staging and port
reporting retain their existing bounded allocations. This one-way receiver
does not replace the original reciprocal radiation load or establish
material, radiation-fit or continuum accuracy.


### Continuous loaded modal and cavity response

`CoupledModalSystem::harmonic_response` solves the original continuous linear
network in peak phasors with `exp(-i omega t)`. Original component damping,
complete signed connection columns, Maxwell material relaxation and attached
multiport radiation all enter the same complex system. Maxwell histories are
eliminated exactly; acoustic oscillators remain explicit, preserving a finite
joint antiresonance at a bare lossless acoustic pole. The existing cancellable
complex LU supplies the solve, with an original-equation backward-residual
gate. Returned original modal and memory coordinates include cycle-mean storage,
source power, each physical loss and their power residual. These are the AC
particular solution about equilibrium, independent of accepted transients and
DC rest lengths; analysis never advances the clock or changes a ledger.
Prepared constant-loss connections reuse the original `g(K)b` damping columns,
including their full cross-mode reaction, while springs keep their original
columns. Returned separate spring/damping efforts and the summed generalized
modal reaction retain those distinct directions. The legacy scalar effort sum
is a channel diagnostic when the columns differ. This is the continuous response
of the admitted prepared image; no new loss law is evaluated at the drive
frequency. G1 tests compare resonant amplitudes with the anchored one-mode law,
and a three-clock ladder converges to the same response through the original
prepared time owner with a nonproportional two-mode connection.

`CoupledModalSystem::harmonic_responses` shares one dynamic-stiffness assembly
and factorization across 1..64 peak-force right-hand sides at a common forcing
frequency. It reuses the original cancellable complex LU solves, original-
equation residual gate, observer and physical accounting for every column.
For d mechanical-plus-air coordinates, G constant-loss links, C connections,
M Maxwell branches and r inputs, the admitted work is d^2*[d+G+r*(C+M+3)];
r=1 retains the scalar allowance. Retained response arrays and struct headers
are bounded at 8 MiB with checked arithmetic before allocation. No coordinate
is truncated and no accepted state changes on success or refusal. Both friction
owners expose the same released batch path. Four production batch tests pass
natively for scalar parity, coherent cross terms, reciprocity, phase cancellation,
work/storage bounds and unchanged state. The complete harmonic module reports
13/14 passing; its unchanged two-body test fails during initial time stepping
before any harmonic call. The same heat-resistance setup failure was reproduced
with pre-batch source c72b13d1881ace6718aefc3c4937ef2a5358b25f.

`ModalPreparation::HarmonicOnly(PreparedModalBudget)` constructs these same
continuous equations without a coupled time propagator, quadratic heat-integral
matrices or held-step friction compliance. `Ordinary` and `Prepared` preserve
the existing playback choices; conversion from `Option<PreparedModalBudget>`
preserves existing callers. The shared preparation owner retains every original
component, connection, material history and radiation pole, and prepares the
same constant-loss `g(K)b` directions using the original stiffness/eigen owners.
Coefficient, material, radiation and internal-component rebuilds preserve the
chosen purpose. Time staging, prediction, compliance and joint-DG time
preparation refuse for harmonic-only images before accepted state changes.
Explicit mode, byte and work allowances still apply; the construction allowance
includes conservative harmonic workspace and coefficient-preparation estimates.
Original physical component initialization and shell reduction still run.

Connection and material rank forms assemble on their exact nonzero support,
using one reused ordered scratch vector. Only exact zero coefficients are
skipped; there is no numerical sparsity threshold or physical truncation. All
retained coordinates, the complex LU, original-equation residual and conservative
work bound are unchanged. Two native support checks compare independent dense
outer products; two native construction checks establish bit-for-bit response
and constant-loss parity with prepared playback, including coefficient updates,
and state-preserving refusal/cancellation. Full-bank and friction construction
regressions remain authored but unrun.

`SharedCavity::harmonic_response` maps original surface forces and incident
neck pressure into that same solve and returns physical internal pressure,
surface volume rate, outward neck flow, mutual exterior reaction and gas
compression availability. Internal and radiated pressures remain reactions,
not extra external drives. Thermal pressure uses the admitted compression law;
viscous and thermal queries refuse outside their sampled admission bands.
Radiation validity remains with its producer, including the declared band of
geometry-derived fits.

A caller work allowance bounds `d^3+d^2*(connections+constant_loss_links+material_arms+3)` before
dense allocation, with `d=mechanical_modes+acoustic_poles <= 1024`. No mode or
pole is truncated. Nonpositive frequency, singular or unrepresentable response,
and missing prepared constant-loss columns explicitly refuse.
This continuous response is not the sampled held-force transfer. Independent
physical pressure/flow condensation, pressure-diffusion availability, separate
loss channels, signed basis invariance and unchanged accepted histories are
tested alongside analytic loaded resonances.

### Continuous finite side-branch resonance

`CoupledModalSystem::harmonic_response_with_load` and the shared-cavity wrapper
retain the main DG bore, every original finite `DgBranchLoad`, and the complete
original MIMO terminal in one continuous frequency solve. Each junction adds an
explicit flow unknown with equal original pressure traces and opposite inflows.
No isolated duct impedance is inverted, including at a bare branch pole.
Terminal channel IDs retain their authored order; a branch channel observes its
actual distal outlet. Original wall memories, material bands, mutual terminal
entries and strictly proper terminal states are preserved. The complete state
and junction count consumes the existing dense work allowance.

Returned branch responses carry original fields, storage, wall losses and both
junction powers. Whole-network accounting sums every owner once. Independent
original pressure-trace checks accompany the complete-equation residual; queries
never advance accepted histories. Analytic cylinder-network peak and linewidth
refinement, a bare branch pole, mutual/strictly proper terminals, branch and
channel permutations, viscothermal bands and shared-cavity accounting are covered
by the native focused tests. These are zero-volume point junctions joining finite
plane-wave ducts, not a resolved three-dimensional entrance-flow model. Finite
gap channels and driven exterior mouths remain explicit refusals.

### Continuous finite side-branch evolution

The explicit exact downstream backend also retains every original branch state
on its common clock. In original DG energy coordinates, junction columns are
`H=[-B_bore,+B_branch]`; the complete normalized Gram solve constructs
`P=I-H(H^T H)^-1 H^T`. The generator `P L P` and projected external columns
preserve the original pressure constraints from the admitted cold state. This
eliminates only zero-storage junction flows; it does not reduce the acoustic or
wall-state basis. The same cached exponential and quadratic-integral owner
evolves the coupled modal, cavity, duct and original rational terminal states.
Terminals admit conditioned invertible D or exactly zero D. For zero D,
original states obey `x'=A x+F U+b u` and `H x=0`; the conditioned solve
`(-H F) U=H A x+H b u` removes only zero-storage flow multipliers. Full
projected inlet coupling, original pressure/work observations and all DG/rational
coordinates are retained. Initial and endpoint constraints and their
differentiated identities must pass; no artificial resistance or state reset
is applied.

Each original duct's unprojected loss form is integrated separately. Junction
flow includes its acoustic-state, modal-velocity and complete terminal-feedback
terms. Both original pressure impulses, both actual integrals of `p(t) Q(t)`,
and swept volume are independently observed. Original per-duct energy closure,
endpoint and integrated pressure continuity, junction work, terminal work and
whole-system budgets must pass before the existing `DgPortTrial` publishes any
history. Terminal work remains signed external exchange because a generic
rational realization supplies no physical storage certificate. Tests compare
24/48/96 kHz constant-drive trajectories, reversed branch ownership, original
wall memory, an independent physical mutual-terminal heat ladder, and late
refusal/cancellation with identical retry. Preparation retains the existing
256-dimensional complete-state allowance; finite gaps, driven mouths,
mixed-rank singular D and unresolved higher-index pressure constraints remain
outside this backend.

### Shared cavity exterior loading

`SharedCavity::apply_modal_impulse` admits simultaneous surface momentum and
exterior neck pressure impulses on the continuing original network.
`apply_pressure_impulses` maps one uniform exterior pressure integral per
original surface through its retained signed volume shapes. Positive pressure
acts inward; the physical neck momentum change is `-area * integral(p dt)`.
Source work is the signed pressure integral paired with mean outward volume
flow. Projection overflow/underflow and unrepresentable resulting flow refuse
before the native momentum, energy and endpoint-force checks publish a jump.
Staging reuses retained scratch and allocates no new trajectory.

Only physical surface and neck velocities jump. Position, cavity volume and
compression pressure, viscous auxiliary velocities, thermal/molecular memory,
proper radiation state, clocks and preceding interval port records remain
unchanged. The native impulse receipt carries the event work separately;
there is no invented heat, force duration or repeated internal pressure drive.
This is the admitted finite-inertia small-signal image, not a shock-front model
or an impulse reset for attached distributed ducts/algebraic constraints.
`tests/shared_cavity_impulse.rs` checks physical mass/work and independent
rank-one cavity motion at two rates, continuing viscous/thermal/radiation
history, signed multimode pressure projection and atomic refusal/retry.

`cavity_modal::shared::SharedCavity::with_exterior_radiation` connects the
existing passive multiport model to the original surface/cavity/neck exact
propagator. A selected neck column is in physical pressure/volume-flow units;
its power-preserving pullback is `Q = S/sqrt(M) q_dot`. Other selected columns
retain the original surface modal velocities, including signed mutual coupling.
The neck must declare `ExteriorMouth`, so exterior mass and resistance are not
counted twice. Additional prescribed exterior pressure remains an incident
source. `neck_exterior_pressure_pa` observes the accepted internal reaction;
feeding that pressure back as an incident source would double-count it.
Analytic coupled-pole motion, exact time partition, energy balance and rejected
step preservation are tested in `cavity_modal::shared::tests`. Geometry/band
validity of the supplied impedance remains with its producer; this composition
does not claim a new BEM fit, distributed neck field or experimental validation.

`SharedCavity::new_with_viscous_neck` replaces the frozen internal neck
resistance with the existing positive `ZkDiffusive` circular-bore law:
Poiseuille drag plus passive `R s/(s+p)` arms, represented by free auxiliary
masses in the same coupled exponential. Only actual bore length contributes
viscous drag and storage; inviscid bore and inner-end inertia remain counted
once, and `ExteriorMouth` radiation can be attached separately. Auxiliary
states consume the original mode, connection and work budgets, retain their
kinetic storage, receive no external port or radiation column, and participate
in accepted-state rollback. `ViscousNeckOptions` declares the sampled Bessel
accuracy band and error budget. This does not claim thermal wall compliance,
distributed neck propagation, nonlinear flow, measured fit accuracy or realtime
throughput. Tests compare independent duct Bessel impedance, a physical-flow
RK4 reference, energy balance and exact timestep partition.

`SharedCavity::with_baffled_exterior` derives the complete exterior reaction
from supplied original signed normal-velocity traces and the actual circular
neck opening in one infinite rigid baffle. The generated annular-sector area
quadrature retains the original mouth radius and area; body traces must retain
every original surface coordinate and reproduce its signed cavity volume
projection. The existing power-consistent Rayleigh owner supplies all mutual
surface/mouth loads, and the existing positive-residue fitter retains its
unchanged training/held-out complex and resistance error gates. Fitting uses a
common mechanical mass-normalized basis before the neck is converted to the
physical pressure/volume-flow convention of the existing attaching method.
Explicit cell, frequency, sweep and complete-pole budgets apply without field
truncation. Preparation is cold; ordinary stepping keeps the existing loaded
propagator and accepted-state accounting. A three-rung mouth quadrature checks
the independent flanged inertance and compact radiation-power limits; signed
mutual loads, energy and exact time partition are also tested. This is a
sampled, spatially discretized infinite-baffle reference, not a curved/unbaffled
body solution, quadrature-error certificate or out-of-band guarantee.

### Cylindrical cavity thermal compression

`SharedCavity::new_with_cylindrical_thermal_wall` uses an explicit physical
cavity radius and the original volume, with an isothermal cylindrical sidewall,
insulated endcaps and the existing uniform-pressure assumption. Positive ZK
thermal compliance is inverted through the existing Jacobi owner into
isothermal equilibrium compression plus positive Maxwell arms, retaining the
adiabatic instantaneous limit. Thermal history, acoustic availability and
diffusion dissipation remain in the original constitutive-memory/exact account.
Pressure and compression-energy observations include accepted thermal branch
reactions once. Relaxed versus instantaneous initial compression is explicit;
viscous-neck and exterior radiation loading compose without changing physical
port count. Three-rung independent Bessel and physical pressure/diffusion
time-domain tests cover this compact reference. No arbitrary-cavity geometry,
finite wall conductivity, distributed acoustics or full wall-temperature
evolution is inferred.

`SharedCavity::new_with_cylindrical_molecular_wall` and
`CylindricalThermalCompression::new_with_molecular_relaxation` bind the
original gas reservoirs to those same radial heat coordinates before lowering
compression. The complete material gas snapshot must match the cavity source.
A positive capacity/conductance pencil yields positive Maxwell arms in one
joint symmetric eigensolve; the existing coupled propagator, harmonic solve,
pressure observation and material-memory work account consume them unchanged.
There is no second heat deposit. Setup explicitly admits the joint coordinate
count (thermal poles times one plus nonzero molecular reservoirs; ceiling 128).
The frozen molecular instantaneous stiffness and equilibrium-gas adiabatic
reference have distinct accessors; the fully relaxed limit stays isothermal.
A 33-frequency check bounds sampled complex and dissipative realization error
to 1e-8 within the original wall band, separately from equilibrium Bessel
admission. Joint spatial error remains a thermal-pole refinement question.
Independent tests cover a quadratic constitutive oracle, three-rung complex
Bessel convergence, reservoir splitting, physical temperature ODE storage/heat,
the original loaded two-coordinate solve, and accepted-time work closure.
This preserves local Debye reservoirs and the original constant Fourier
conductivity; separate molecular-energy diffusion, finite wall conductivity,
nonuniform pressure, flow and changing background state are outside this lane.

### Molecular cavity compression

`SharedCavity::new_with_molecular_relaxation` maps the original material
`RelaxationState::bulk_relaxation` spectrum through the actual cavity volume
and neck area into the existing reciprocal modal solve. Internal pressure,
compression storage, molecular loss, harmonic response and accepted histories
use the same source gas. Original surfaces, neck viscous memory, downstream DG
and exterior radiation retain their owners and work accounts. The supplied
material law must match the complete gas snapshot; missing/extra history,
setup-budget failures and cancellation refuse publication. Independent
internal-temperature dynamics and thermodynamic availability check the
time-domain pressure/energy; direct caloric response checks harmonic pressure
and dissipation in `cavity_modal::shared::molecular::tests`.

The medium is frozen, the compact pressure spatially uniform, and the cavity
wall adiabatic. Molecular loss is acoustic availability dissipation, without a
second gas-heat deposit. For a conducting cylindrical wall, use
`SharedCavity::new_with_cylindrical_molecular_wall` above: its joint caloric law
replaces the adiabatic specialization. The two compression owners cannot be
added as independent shunts. No measured instrument fidelity or evolving gas
chemistry is inferred.

### Accelerating material-fiber supports

`MultiBowPerformance::with_fiber_support_program` supplies declared continuous
piecewise-linear axial velocity to the existing `with_fiber_support_motion`
owner. That owner projects `Gamma_n = integral(rho A phi_n dx)`, supplies
`-Gamma_n a` to the original modal/contact/nonlinear solve, uses step-average
support speed, and migrates the crossing along the original material hair.
There is one inertia law and one accepted support state. The program rejects
competing tangential gesture speeds or another support-motion control. Velocity,
acceleration and complete translation extrema are bounded before playback;
the original finite span independently bounds the moving crossing. Original
materials, pressure, state and sample clock survive; callback partition cannot
change the motion.

The native mechanics stores fiber-relative kinetic energy. `FiberSupportFrame`
adds `V Gamma.qdot + M V^2/2` to recover laboratory mechanical energy. Physical
support work is native bow drive work plus inertial-coordinate work plus the
change of that kinetic correction. Both direct and programmed support motion
use this account. Component-state and work bounds admit representability before
any physical step; publishing an accepted frame is infallible. Inertial work is
never heat. Prescribed normal load and dry friction remain this moving-hair
rung's assumptions; normal/rosin/thermal history transport and stick flexure
are not supplied. Focused regressions compare free acceleration with the
independent sine-wave solution, work closure, width refinement, material
substitution, callback partitions and the original moving-contact owner.

### Coupled extension and normal deflection of material fibers

`TensionedFiberRibbon::with_axial_transverse_strain` uses the original
`fs_nlmodal::longitudinal::LongitudinalString` storage, with each quadrature
strip's original mass, EA and tensile preload. Normal deflection changes
`u_x + y_x²/2`, generating axial waves and changing local tension and contact
reaction. All original fiber coordinates and linear material losses remain in
the shared exact-linear/contact/discrete-gradient owner. Its existing axial
quadratic storage is subtracted once from the nonlinear correction. Disjoint
cold material potentials compose under one identical solve configuration;
overlapping coordinates or incompatible budgets refuse. Existing string, body,
cavity and pressure states are retained. Continuous-span slope/strain bounds
and positive local tension gate every trial. This fixed-end moderate-slope
rung retains selected axial inertia; omitted axial harmonics equilibrate
statically. Axial/normal modes and time require refinement. It adds no stick
camber, moving contact, nonlinear viscous cross terms or real-instrument claim.
Focused G1/G3 checks compare independent continuum strain energy and gradient,
quadrature-width invariance, normal-to-axial energy transfer, and rollback.

### Hereditary response of material fibers

`TensionedFiberRibbon::with_material_relaxation` binds the original section's
layerwise generalized-Maxwell Young spectra to both axial and bending motion.
Each source arm projects to `DeltaE*A*k^2/mu` and `DeltaE*I*k^4/mu`, respectively;
each width cell retains its own complete first-order histories in the existing
coupled string/body/contact owner. Material times and layer identities are not
averaged, and no fitted modal Q replaces the hereditary strain. Width refinement
partitions fiber count and mass together. The exact same source section,
complete coupled instantaneous band, physical clock and history/setup budgets
are admitted before playback. Existing body/string material history survives
the appended fiber channels, and support-motion/contact rebuilds retain them.

The source law remains fixed-state linear Maxwell viscoelasticity. With
axial/transverse strain selected, axial arms act on the same retained
`u_x+y_x^2/2` through `fs_nlmodal::longitudinal::relaxation::StrainRelaxation`.
The shared material owner integrates each endpoint-strain ramp exactly;
interval-mean stress is pulled back at midpoint mechanical coordinates and
physical heat is integrated independently. This keeps nonlinear geometric
work and original source memory together. Linear axial and bending histories
join the same exact augmented mechanical propagator and also admit fast decay
without a midpoint material-time restriction. Source-selected Kelvin--Voigt viscosities
remain explicit independent parallel losses. Nonlinear mechanical motion,
material bandwidth and strain require their original resolution checks;
exact constitutive integration does not imply exact coupled trajectories.
Preload creep, thermal evolution, stick flexure and measured horsehair
parameters are not supplied by this binding.

The full-strain longitudinal string binding uses this same constitutive
update and quadratic-strain work pullback. Its source admission no longer
requires dt/tau <= 2 for these exact ramp histories; the complete instantaneous
mechanical frequency, original material support, strain bounds and all other
participants' clock limits remain binding. Linear axial branch selection uses
the exact coupled exponential with these same source and mechanical limits.
Synthetic layered-continuum tests compare complex displacement and dissipated
power with independent Young-modulus assembly, and exercise source substitution
through real contact, support reversal, accepted work and callback invariance.

### Material gas loading in the file-driven plate renderer

`observer rayleigh X_M Y_M Z_M REFINEMENT_LEVELS MAX_HZ` selects distributed
surface observation in native plate performances. One original eigensolve
retains the complete signed P1 displacement trace alongside its unchanged
modal mass, damping and nodal force projection. Positive triangle subdivision
feeds `pcm_wav::baffled::BaffledPressure`; all original modal velocities,
including zero-net-area modes, reach its existing bandlimited derivative and
path-dependent Rayleigh kernels. No new wave kernel or frequency retuning is
introduced. The fixed baffle-frame microphone, declared gas, acoustic phase
and clearance controls, bounded surface count and retained natural-frequency
band are admitted before receiver allocation. Cold preparation and bounded
convolution preserve callback partition independence. Ensembles align this
receiver's extra numerical latency while retaining the original physical paths.
Physical path delays and numerical filter latency are reported separately; no compact r/c delay
or compact atmospheric filter is added to this observation. Existing optional
piston reaction remains an independent mechanical choice. The approximation
is a stationary P1 surface in an infinite baffle with lossless exterior gas;
source bandwidth, structural/modal convergence, acoustic quadrature, absent
scattering and noncompact loading retain their own physical limits.

The optional `kinetic-mass lumped|consistent-dkt` record after `mechanics`
selects the actual assembled kinetic operator. Omission retains legacy
arithmetic. `consistent-dkt` uses the exact P1 transverse integral and the
original quadratic DKT rotation-field integral from `fs-plate`, with each
triangle's authored or source-resolved density and thickness. The existing P1
force and signed-area pressure observations retain that same transverse field.
`certified_chart_radiators_with_mass` exposes this choice to code consumers;
modal frequencies, normalization, force response and any piston reaction all
use the resulting mass matrix. Integration occurs only during preparation,
with no added callback cost. The input hash, `PlatePerformanceInfo` and CLI
sidecar retain the choice. Native file/WAV regressions check changed physical
output, legacy parity and block-independent replay. This remains thin-plate
rotary inertia, without transverse shear or a spatial-convergence guarantee.

The optional `radiation-load none|baffled-piston` record after `observer` in a
plate performance selects the existing compact piston reaction. The one
resolved ambient gas and output clock supply both mechanical loading and the
retarded pressure observer. Physical face area, modal mass and signed volume
and force projections remain those of the original chart reduction. The input
hash and CLI sidecar retain the load selection. Omission preserves the original
unloaded image; no material loss or radiation model is inferred.

`CompactPlateVoice::from_radiators_with_piston_loads` exposes the same path to
code consumers. A retained filter owns its original gas snapshot and sample
clock: mismatched consumers or step durations refuse before filter mutation;
reattachment of the identical load preserves its history. Gas snapshots pass
the material owner's coherence admission. This remains the existing sampled,
per-mode compact piston approximation, without mutual radiation, atmospheric
absorption, a whole-fluid energy claim or continuum/experimental validation.
G3 regressions distinguish a pressure-induced change in mechanical response
from density-only output scaling, check the actual retained body's kinetic
energy, and compare loaded file/WAV renders across callback partitions. The original piston
decay/frequency oracle retains its tolerance and uses a coherent low-pressure
gas for its vanishing-load limit.

### Material-driven finite-body sliding

`tribothermal::dynamics::SlidingPair` adds force-driven tangent inertia to the
original material sliding contact and finite enthalpy network. Each complete
uniform body's original thermal-cell mass also supplies its inertial mass.
Both velocities and applied forces have two components in the original fixed
tangent frame, ordered by body A then body B. Centre-of-mass momentum and two
Newton--Krylov relative-velocity residuals determine the complete planar motion;
every trial restarts the original contact/thermal candidate,
including its actual endpoint temperature, pressure, speed and slip-history
source queries. Separate impulse and whole-system energy gates precede joint
publication of velocities, original slip/work history and finite temperatures.
Original friction work is an internal transfer; external force work is measured
as the dot product with each body's endpoint displacement. Physical friction alone is partitioned
as heat, while backward-Euler kinetic loss remains an explicit numerical account.
The fixed patch/normal-load reduction passes released motion and inertial
compliance `dt²(1/mA + 1/mB)` to the original tangential return map. That owner
selects exact sticking or vector sliding; the material source queries the
Euclidean speed, and accumulated distance follows actual slip arc length.
The source motion converges to that same two-component native slip. A stuck
branch queries exact zero speed
and leaves slip/work history unchanged. Bodies can accelerate together, stop
and restart in the opposite direction, under the native endpoint static-first
convention. It does not locate the physical stopping instant within a coarse
interval: a backward-Euler stop may have numerical kinetic loss and zero endpoint
friction heat. Refine time to resolve stopping heat; numerical loss is never
thermalized. Zero-speed material support is required. No normal, rotational or
contact-geometry evolution is added. Focused synthetic regressions cover static
capacity, signed source feedback, oblique motion, rotation of the tangent frame,
stopping/restart, continuation and rollback;
native execution remains required for qualification. The internal same-trial
friction callback allows other energy owners to share this original solve;
compliant consumers must gate queried motion against the native solved slip.

The `SlidingPair::from_material_bodies` constructor also binds the original
`fs_mbd_material::ResolvedBody` objects to their ordered interface materials
and finite caloric cells. It requires uniform same-card bodies with matching
mass and total reference volume. Motion and kinetic storage consume the
retained canonical body inertias; accepted candidates share those immutable
sources. A density change at fixed geometry therefore changes both inertial
response and the matching finite heat inventory through data resolution.
Mismatched identities, masses or geometry refuse before motion advances.

`SlidingPair::step_planar_candidate` advances turning and stopping motion in
that same fixed tangent frame through the original material/contact/enthalpy
candidates. Each segment uses midpoint displacement, with a scalar radial
speed solve for the current isotropic friction law. Physical traction work
heats the original finite bodies; no backward-Euler kinetic defect is booked
as heat. Static reaction is considered only at an accepted zero relative
velocity. A stopping-time root must also satisfy the full vector impulse
gate; unresolved turning is subdivided instead of accepting a mirror-velocity
midpoint root. Both bodies' accumulated impulse defects and the whole-call
energy account must pass before any segment is published.

Constant collinear Coulomb motion is exact at arbitrary tangent orientation,
including arrest and restart. General curved paths require refinement of the
declared maximum segment duration, independently of root accuracy. An oblique
limiting stop is approximate within its vector impulse/time budgets. G1/G3
tests compare curved-stop integrals and a separate continuous RK4 turning
reference on three time rungs, and exercise rotated frames, distinct static
and kinetic capacity, source temperature feedback, and transactional refusal.
Temperature-dependent coefficients still use endpoint thermal states and
segment mean speed; thermal transport remains first order. No normal motion,
anisotropic material law, within-stick thermal event localization, or real
material calibration is inferred.

### Directional material friction with finite heat and inertia

`tribothermal::directional` retains the original four directional coefficient
claims, ordered material/phase identities and finite caloric cells. Each thermal
trial refreshes declared surface/mean temperatures, nominal pressure, sliding
speed and unsigned path; undeclared source conditions stay fixed. The native
associated ellipse supplies actual traction/slip work, partitioned explicitly
between the two bodies and any named outside destination.
`DirectionalPair::from_material_bodies` also retains the original source mass,
density and geometry receipts.

Force-driven planar candidates use midpoint compliance `h²/(2 m_reduced)` and
the native scalar ellipse projection. Only speed/path-dependent sources need
an additional bounded scalar consistency root. Accepted source path never
exceeds actual accumulated motion. Static capacity applies only from accepted
rest; moving segments use kinetic coefficients. Full-vector arrest, turning
subdivision and rank-one partial arrest retain motion along a frictionless
material axis, including widely separated component scales. Independent
per-body impulse and whole-system work/heat gates precede joint publication;
cancellation or exhausted budgets leave all original owners unchanged.

Fixed normal load, nominal area and tangent/material frame are required.
General turning and temperature/speed/path dependence require timestep
refinement; endpoint thermal transport is first order, and arbitrary weakening
is not claimed to have a unique solution. No normal, rotational or
evolving-patch dynamics or experimental calibration is inferred. Nine focused
synthetic test groups cover source ownership, finite inertia/heat, independent
continuous turning refinement, static/kinetic separation, temperature and motion
feedback, partial arrest, monotone continuation and atomic refusal.

### Material contact loads on a field actuator

`tribothermal::field_motion` closes the existing dry-contact active-set/enthalpy
owner with `MaterialFieldMotion`, including optional sourced axial support.
Two displacement residuals identify actual field gap change, original solved
contact slip, and the motion used to query interface material data. The native
traction supplies the opposing field force. Field/circuit/body/support history
and contact/thermal history publish together; inertia is not integrated twice.

The same thermal trial resolves bulk electrical transport, pinned electrical
contacts, and pinned mechanical friction. Original conductor/contact losses and
the original friction heat partition enter explicit receivers exactly once.
External thermal supply is independently evaluated before either internal loss.
Field, kinetic, support elastic, interface and thermal storage plus exported
friction work must close original electrical, mechanical and thermal supplies.
A single descriptor budget includes all discarded outer and thermal trials.

Surface A must match the single-region field body's source card, mass and
reference volume; surface B is a stationary finite thermal guide. Normal force,
area and tangent frame are fixed. Static reaction at rest uses an exact zero
speed/history source query; uninterrupted sliding in either direction is
supported. An unresolved stop/reversal refuses;
there is no velocity clamp, impact impulse, moving thermal mesh, or new normal
contact law. Focused `tribothermal::field_motion` regressions cover independent
constant-charge trajectories, both static reaction directions, body/source
inventory matching, authored/callback external heat, caloric/resistance/friction
feedback including electrical contacts, cancellation and global-budget refusal.
Syntax parsing and independent numeric checks passed; native execution is pending.

### Typed gesture consumer bindings

`gesture_path::bind_paths` binds each description-validated parameter path to
an explicit existing pressure voice or actual material-string bow member.
Every source path requires a destination; incompatible/unconsumed paths return
`GestureError::NoConsumer { path }`. `BoundPaths` retains the original resolved
v2 schedule, alias-migration identity, description and exact destination map for
render provenance. V2 schedules remain typed throughout lowering and call the
original pressure/friction sampling loops without changing time rounding, ramps, station
checks, integration, cancellation, force transfer or stored physical state.
Numeric u8 alias decoding is isolated in fs-scenario's compatibility module;
typed destination members use usize without a legacy u8 conversion.
The existing single-string, material-bank and simultaneous-bow entry points have
path-bound counterparts returning that provenance alongside their renderer.
Legacy file entry points remain v1 compatibility adapters; a typed schedule
requires description admission and explicit path bindings. No new consumer or
silent fallback is invented for topology events, frets, valves, pitch, prestress
or terminations. G0 `gesture_schedule::unconsumed_targets_reach_a_typed_consumer_or_refusal`
checks every alias and bitwise pressure-control equivalence; original render
equivalence tests remain unchanged. Addressing alone establishes no physical
validation or conservative topology-change claim.

### Declared per-note control resolution

`description::per_note` resolves `MinimumThresholdPressure` using a model-only
monotone onset predicate and explicit SI domain, pressure tolerance and maximum
evaluations. Its API takes no audio/reference input. The non-oscillating lower
and oscillating upper bracket, every probe, model identity, objective, target,
note, domain and budget enter a deterministic content-addressed receipt.
`record_resolution` binds the receipt to the corresponding Performance control;
`fs_scenario::performance::lower` retains it in render provenance. Other declared
objectives remain representable but refuse this pressure-specific resolver.
Model identities must bind initial conditions and the predeclared observable;
evaluation cancellation is the model's responsibility. The search is bounded
and synchronous and publishes only a completed result. Monotonicity is an
explicit model obligation; a sampled bracket is not a proof of global
monotonicity or asymptotic bifurcation onset. G1 `per_note_rule::pn_001` uses
the existing reed fixture's 0.25 s sustained-pressure observable with 18.75 Pa
bisection width; G0 `pn_002`–`pn_004` cover API shape, provenance and refusals.
The authored four-note `music_lane.sh --fixture midi-reed` exercises MIDI,
lowering, model-only resolution, schedule rendering, replay and named refusal.
This fixture is plumbing evidence and makes no reed realism claim.

### Continuous force ramps through the complete modal network

`CoupledModalSystem::into_linear_forcing` prepares fixed-clock continuous linear
force records for every original retained coordinate. `LinearForceModalSystem`
keeps the original component damping, signed spring and dashpot columns,
prepared constant-loss damping columns, Maxwell internal strains and reciprocal
radiation states in one generator. The wrapper retains all accepted history and
provides read-only system access; changing coefficients or attachments requires
consuming it and preparing the changed system again.

External work integrates the actual linear force times velocity, independently
of the energy residual. Each component, physical connection, material and
radiation loss integrates its original positive source squares along that same
trajectory through `fs_phs::linear_exact::PreparedSourceLoss`. Each connection's
port work is minus its spring-energy change and its independently integrated
loss. No endpoint quadrature, total-energy remainder, diagonal modal reduction
or numerical dissipation supplies physical heat.

The original exact block, `LinearSplit`, `PreparedQuadraticFlow`, candidate
admission and publication owners remain authoritative. Additional impulse-valued
force coordinates carry the affine drive without changing physical storage.
Original state, reaction, pressure, energy, endpoint material-domain and
cancellation gates run before the complete history and port account publish.
Held-only material-trajectory queries remain unavailable for a ramp; no finite
thermal receiver or nonlinear held-contact reaction is implicitly reinterpreted.

Cold preparation admits total cache/workspace bytes and multiply terms before
dense assembly. Twice the complete physical-plus-force-record dimension must
fit the shared 256-coordinate positive-integral capacity; larger configurations
refuse without dropping states. Warm steps reuse flow and positive heat factors,
while retaining existing candidate and port-account allocations. There is no
new allocation-free whole-system or measured real-time claim. Focused tests
compare complete motion, four heat channels and signed work with independently
integrated continuous equations, force-ramp subdivision, equal-endpoint held
forcing, free/critical/overdamped scalar physics and transactional refusal.

### Port-schema adoption contract

`PortedModule` exposes stable local schemas, coordinate charts, internal N-way
junctions and the last accepted physical-time transfers for CoupledModalSystem,
SharedCavity, GasPiston, DoubleActingGasPiston, DynamicAperture, ModalAperture and
the ReedBoreVoice characteristic inlet. All schemas pass PortSchema::try_new.
Graph embedding must namespace local identities. The chart timestamp names the
clock origin; each audit binds its own exact accepted interval. Held effort work
uses the original conjugate displacement increment, never logical ticks as time.

Adapters observe the existing solvers; they do not feed back into their arithmetic
or acceptance. Typed transfers and window evidence are materialized lazily on
request. Exact-network connection work is minus the spring-energy change and
the independently integrated damper loss along the same full trajectory, never
endpoint reaction times displacement. Multiple damped columns require bounded
per-column integration scratch; this is not an allocation-free observer.
Rejected trials preserve the preceding account. An independent observation
refusal is retained as `ModuleObservationFailed`, and observational overflow as
a typed transfer error; neither publishes an invented zero record or rewinds
an already accepted step.

SharedCavity retains its mass-normalized modal chart and signed volume projection.
The canonical schema is the explicitly transformed reference-coordinate pair:
F = f_q sqrt(M_ref), v = q_dot/sqrt(M_ref), with M_ref=1 kg, not a specimen modal
mass or a point-force measurement. Raw owner units remain N/sqrt(kg) and
m sqrt(kg)/s. Its compression junction uses the original signed volume pullbacks
to pressure/volume-flow coordinates, including neck and compression storage.
Each CoupledModalSystem column similarly names its B_i^T q_i work charts and
one connection-storage endpoint. DoubleActingGasPiston's two gas faces share a
single common-flow wall junction; internal face work is not an exterior supply.

Available energy windows use the owner's actual storage change, prescribed work,
loss and existing joule tolerance. No balance is reconstructed as a leftover.
Piston heat ports use the source-backed discrete gas-midpoint temperature and
Q/T_midpoint entropy quadrature; gas entropy change remains independently
measured. This energy projection claims neither gas species/element balances,
exact isentropy nor external-reservoir thermodynamic validation. Negative
numerical entropy generation remains a visible WindowAuditReport violation.

MagentaBadger's ruling in Agent Mail 2636 explicitly leaves two audits unavailable:
DynamicAperture has no owner-declared energy tolerance, and ReedBoreVoice /
DelayedFilter has no physical line storage/loss owner or joule tolerance. Their
real transfers remain available, but ModuleWindowAudit::Unavailable names the
missing inputs. Those gaps keep IRP .162 open; unavailable is never a passing
window. No stepping or music-lane golden is intentionally changed.

G0/G3 tests in tests/port_schema_adoption.rs cover admitted schemas for all seven
owners, available/unavailable accounts, sign-flip refusal, two-port parity,
three-component/shared-cavity N-way descriptors, the common-wall descriptor,
and accepted-record preservation on refusal. Each module emits a native-run
JSONL row with residual/tolerance or explicit null/unavailable fields. Native
verification and reed/string baseline hashes remain central-run requirements.

### Integrator doctrine and current audit boundary (IRP M5.2-W1)

Render-path linear parts must use exact held-input propagation or explicitly
frequency-matched energy-consistent propagation. Nonlinear storage retains its
discrete-gradient owner; generalized-alpha is not admitted on render paths.
Held effort exchanges use the same displacement work at both ports, rate
changes require antialias filtering, and phase/contact/solver budgets must be
admitted explicitly. A midpoint path that remains is admitted only inside the
fs-phs one-cent band (`admit_midpoint_band`, pole-rate bound times dt <= 0.04);
the table below says which paths are exact and which are banded.

`tests/integrator_warping_audit.rs` contains the current path inventory and G1
one-step-map measurements at 1/5/10 kHz and 48/96/192 kHz. The common fs-phs
linear owner is measured for its callers; that is not a complete nonlinear or
connected render. Public ImpactSystem and Dynamic/Paired/ModalAperture seams
are exercised directly, including the uniform scalar condensation. The exact
modal owner is the control. Midpoint phase is compared to
`2*atan(omega*dt/2)/dt` within 1e-6 cents. The closure test
`every_render_path_is_listed_with_integrator_and_measured_warping` accepts
a row only if it is exact, admits no more than the one-cent band, is the
first-order network row with its measured order ladder, or is a PianoAudio
factor whose companion measurement stays within 1 cent. A path outside its
band refuses with a typed `TimeResolution` error; nothing is retuned.

| Current path | Linear owner | Current band / audit boundary |
| --- | --- | --- |
| Fixed-end KC strings (`realize_kc_string`, `realize_coupled_kc`, `realize_coupled_ode`) | Prepared exact quadratic propagation with the complete nonlinear discrete-gradient endpoint, admitted using the actual propagation and heat matrices | Retained bare/bending modes fit the observer Nyquist guard; nonlinear/contact motion requires time refinement. Complete models above the existing matrix ceiling refuse without mode truncation |
| Nonlinear plate (`thin_plate::VkBody`) through fs-phs `step_render` | Exact ZOH when the storage is exactly quadratic, otherwise Gonzalez midpoint | Midpoint admits only the one-cent band on the split's pole-rate bound and refuses with `TimeResolution` |
| Dirac joins (`realize_dirac_join` / `step_descriptor`) | Gonzalez midpoint | One-cent band on the descriptor pole-rate bound; typed `TimeResolution` refusal |
| Glottis fold island, coupled ODE chains, blown and reed ODEs, brass lip held-force island, plate-cavity transformer | Exact ZOH (fs-phs `linear_exact`) | No warping; nonlinear forces are held over the step |
| PlateBank complete body/material/cavity/exterior state | Prepared full exact ZOH for quadratic storage, otherwise Gonzalez midpoint | Reciprocal loads share one accepted solve; the original nonlinear pole-rate band and complete-state capacity apply |
| TsDriver speaker coil/cone/panel | Full coupled exact ZOH | Held voltage, independently integrated current/work/loss; dense reference step, not an allocation-free callback |
| ImpactSystem | Gonzalez midpoint | One-cent band on the complete retained block (modal stiffness, volume-spring dyads and resistance), typed `ImpactError::TimeResolution`; the per-mode 0.9*pi guard remains a separate validity check |
| ModalAperture, uniform condensation, and Dynamic/PairedAperture with `BladePropagator::Midpoint` | Midpoint response | One-cent band (pole-rate bound times dt <= 0.04, DEF-04) refuses at construction and at every step; the material-modal adapter additionally caps its caller's bound at 1 |
| PairedAperture with `BladePropagator::ExactZoh` | Exact ZOH per elastic blade (scaled exponential; Van Loan viscous loss) | No warping at any omega*dt (audit rows ~1e-11 cents up to omega*dt 1.31). Contact/Bernoulli forces are held over the step and solved at the midpoint opening; blade material memory refuses |
| DynamicAperture with `BladePropagator::ExactZoh`, authored or plate-bound (closure, Maxwell arms), and every plate-valve performance file | Exact ZOH of the opening plus its Maxwell arms (scaled exponential; Van Loan viscous and arm loss), through `reed_bore::contact::step_with_opening` | No warping (audit rows <= 3e-10 cents up to omega*dt 1.31); arms checked against an independent RK4 integration. Contact/Bernoulli forces held over the step and solved at the midpoint opening |
| CoupledModalSystem | Exact ZOH of the full retained bilateral pencil (component modes, connection springs and dashpots), split into connected K/C islands (cd7f65bc7); nonlinear contacts keep held discrete-gradient reactions around the same affine response | Coupled algebraic Nyquist bound; measured against the exact two-body eigenproblem; the two-body sticking traction matches an independent exact-flow integration within 1e-10 N |
| ModalAcousticTimeModel | Exact ZOH | Caller-declared Nyquist guard; no frequency-warping claim for its observer transfer |
| Source-backed GasPiston | Midpoint gas/momentum | One-cent band on the gas-spring rate sqrt(rho c^2 A^2/(m V)) at the start, midpoint and end states; typed `PistonError::TimeResolution` |
| Circuit/thermal/network piston | First order in general | Observed order 1 +- 0.1 on the four-rung companion ladder (`FS_IRP_NETWORK_AUDIT_LOG`); no exact-linear claim |
| Live PianoAudio at 1/2/4x | Exact local blocks plus finite-step coupled solves | Companion-measured coupled-board frequency within 1 cent, window discrepancy included, at each factor (`FS_IRP_PIANO_AUDIT_LOG`) |

Fixed-end KC members carry every admitted Prony bending memory coordinate in
their exact quadratic block. Their material amplitude decay therefore does not
inherit the former midpoint `dt/tau <= 0.1` or `dt*omega_instant <= 0.1` limits.
The original material frequency domain, relaxed and instantaneous observer
Nyquist checks, complete-state resource bounds and nonlinear/contact accuracy
requirements still apply. Moving-end joins retain their reference clock limits
and their independent coupled-rate admission. G1 consumer regressions compare
5/10 kHz modes at 48 kHz with `dt/tau = 16` against independently refined Maxwell
ODEs, including amplitudes, material history, positive heat and held-force work;
a public assembly regression carries the admitted decay through to pressure.

The coupled-mode test uses the exact symmetric stiffness eigenproblem and a
three-rate phase-error ladder. The coupled pencil steps exactly, so the scalar
midpoint atan oracle does not apply to it. Piston tests reuse the actual
source-backed thermochemical and thermal/electrical fixtures. The speaker
migration has a held-voltage time-partition regression and a public-state pole
audit. The audit and its companions pass natively (M5.2-W1 closed on ts1 b46 at
a4b6c4c22); consumer fixtures that request clocks outside a midpoint band
refuse with TimeResolution and must move to a resolved clock or an exact block.

- `PORT_SCHEMA_VERSION = 2`; `StableId` admits a canonical transport-safe ID.
- `PortValueShape` distinguishes scalar, non-empty vector/tensor, and field
  values; field shapes carry separate neutral `fs-iface::SpaceType` roles for
  effort and flow. `PowerPairing` distinguishes scalar product, Euclidean dot
  product, and field duality with explicit integration-measure dimensions plus
  `FieldMeasureSide` (effort-density or flow-density).
- `PortKind` covers mechanical translation, rotational torque/angular velocity,
  fluid pressure/volume flow, thermal temperature/entropy flow, electrical
  voltage/current, magnetic mmf/flux rate, and chemical electrochemical
  potential/amount flow. Each kind owns canonical generalized effort/flow
  dimensions; semantic kind identity remains distinct even when dimensions
  coincide (for example torque and energy).
- `CoordinateBinding` makes basis, frame, and positive orientation explicit;
  `PortTimestamp` carries a named logical clock and tick;
  `ConservationRole` is canonicalized into a sorted, duplicate-free list.
- `PortSchema::try_new` binds stable ID, legacy/current `PortKind`, six-base
  effort/flow dimensions, value/field shape, coordinates, power pairing,
  timestamp, and conservation roles. V2 admission checks that the contraction
  matches the shape, effort × flow has watt dimensions without exponent
  overflow, the measure-adjusted generalized dimensions match the declared
  `PortKind`, and PR-2 schema-only kinds name their non-energy conserved flow
  where one exists. The original three seed schemas retain their PR-1
  Energy-only role vectors exactly.
- `PortKind::scalar_seed_schema` constructs one canonical scalar coordinate of
  any kind without inventing identity, frame/basis, or clock data; retained
  goldens use it to migrate the existing mechanical, fluid, and thermal seeds.
- `STREAM_PORT_VERSION = 1`; `StreamPort` is not a scalar `PortKind` or a fifth
  thermodynamic relation. It bundles signed mass (`kg/s`), canonically ordered
  species/element amount (`mol/s`), three momentum-rate (`N`), energy (`W`),
  and entropy-rate (`W/K`) values under one `StreamChartBinding`. Its fixed
  roles are Energy, Mass, Amount, LinearMomentum, and Entropy.
- `StreamChartBinding` separates the constituent-basis artifact and its
  explicit, canonically ordered species/element axis from the spatial
  basis/frame/orientation. It also binds the state schema, chemical reference
  state, logical clock/tick, gravity datum, and the closed
  `StreamStressWorkConvention`. Stream admission is owner-outward only until a
  public pullback exists.
- `StreamEnergyChart` structurally selects exactly one of: the canonical
  moving-stream enthalpy chart
  `mdot * (h + |u|^2/2 + g*z) + W_deviatoric`; internal energy plus pressure and
  deviatoric Cauchy work; or one coordinate from an exact mixture
  Euler/Legendre family. Internal-energy, enthalpy, Helmholtz, and Gibbs
  selections all reconstruct canonical enthalpy from the retained conjugate
  terms before transported energy is formed. The caller-declared stream energy
  rate must equal the selected chart bit-for-bit.
- `PressureWorkCrosswalk` recomputes `h = e + p/rho`,
  `volume_flow = mdot/rho`, and `mdot*(h-e) = p*volume_flow` exactly.
  `EulerLegendreCrosswalk` recomputes the mixture Euler identity plus enthalpy,
  Helmholtz, and Gibbs transforms. Both also require an
  `ExactIdentityProofRef` bound to the complete stream context.
- `ChemicalEnergyAccounting` admits either energy embedded in the selected
  state potential or one separately proved species-potential power term. It
  refuses dual ownership, a foreign chemical reference, and a foreign or
  wrong-kind partition proof. Explicit species-potential and Euler/Legendre
  modes require a species-only constituent axis.
- `ConservativeJunction::iterate_added_mass_fixed` /
  `iterate_added_mass_aitken` are schema-bound migration bridges for the
  retained mechanical scalar fixture; they must remain bitwise equal to the
  legacy results and are not general FSI operators.
- `ConservativeJunction` admits two conjugate schemas and evaluates the finite
  scalar seed as `SchemaInterconnection`/`SchemaPort` (shared effort, opposite
  flow). It localizes metadata mismatches and refuses non-scalar/non-finite
  evaluation.
- Four closed primitive variants are distinct in `PortPrimitive`:
  `ConservativeJunction`; `StorageElement` with Hamiltonian/free-energy state
  schema and constitutive-gradient operator; `DissipativeRelation` with a
  constitutive family plus mandatory monotonicity/nonnegative-production
  evidence reference; and `SourceOrReservoir` with an explicit signed
  `AccountingBoundary` carrying the same basis/frame/orientation as its port
  plus included-source/external-reservoir treatment.
- `ACCOUNTING_WINDOW_VERSION = 1`; `AccountingWindowInterval` names an exact
  nonempty interval in one logical clock without treating tick differences as
  seconds. `WindowEvidenceRef` binds a receipt/verifier/statement digest to
  that interval and to one exact `WindowEvidenceRole`: manifest closure,
  manifest row and local port, initial/final inventory, integrated row and
  local port, boundary temperature, element projection, or charge projection.
- `WindowElementSchema` either declares a nonchemical no-claim or retains a
  nonempty canonical element axis plus interval-bound projection evidence.
  Stream and chemical-power manifests require the audited form.
  `WindowChargeSchema` requires direct-coulomb, species-to-coulomb, or
  proven-neutral evidence; neutrality admits only exact zero charge values.
- `WindowManifestEntry` retains contribution, exchange, local port,
  counterparty, topology role, port kind, logical timestamp, optional external
  `AccountingBoundary`, and row-binding evidence. Stream entries additionally
  retain the complete `StreamChartBinding`; every stream in a window must use
  the same state schema, constituent basis/axis, chemical reference, spatial
  coordinates, logical clock, gravity datum, and stress convention. Internal
  exchanges are admitted only as two reciprocal rows of one kind and timestamp.
  Manifest-closure evidence carries the complete canonical row subjects,
  including topology, boundary, and retained stream context, so it cannot be
  rebound to a different or truncated same-interval manifest.
- `WindowBalance` and `WindowInventorySnapshot` carry finite typed energy,
  mass, canonical elemental amounts, charge, and entropy at an exact endpoint.
  `IntegratedWindowTransfer` uses one convention for every signed value:
  positive into the audited system. `BoundaryEntropyBreakdown` makes advected,
  diffusive/chemical, heat, and radiation terms mandatory and retains a
  contribution-bound constant/profile temperature reference whenever heat or
  radiation is nonzero.
- `AccountingWindowSpec` owns the exact manifest, a common stored/stream energy
  reference contract, element/charge/entropy conventions, typed nonnegative
  tolerances, optional reference environment, and closure evidence. An empty
  manifest is legal only as an explicitly evidenced isolated window. Filled
  rows must cover the declared manifest exactly; ordinary power rows may carry
  only balance axes compatible with their `PortKind`, while stream rows carry
  the full multi-balance vector.
- `WindowAuditReport::audit` returns `Err` only for malformed, incomplete,
  foreign, nonfinite, or semantically incompatible input. A physically failed
  audit returns `Ok(report)` with retained `WindowAuditViolation`s. It computes
  `R = (final - initial) - sum(transfer_into)` for energy, mass, each element,
  and charge; computes `S_gen = delta(S) - sum(S_boundary,into)`; applies
  absolute typed tolerances to equalities and a one-sided lower bound to
  entropy generation; and retains canonically ordered per-port and per-exchange
  ledgers. Every internal exchange is independently gated for zero net energy,
  mass, each element, charge, and entropy transfer, so defects cannot cancel
  between pairs.
- `ExergyEnvironment` and `ExergyLedger` optionally report the narrow
  Gouy-Stodola diagnostic `T0 * S_gen` for a finite positive explicit `T0`.
- The raw scalar `Port` remains a backwards-compatible, non-admitting numeric
  container. The legacy `conjugate_to`/`interconnect` migration oracle composes
  only the original mechanical, fluid, and thermal seed kinds and refuses all
  schema-only kinds; raw construction or arithmetic alone carries no
  dimensional, coordinate, clock, identity, or conservation certificate.
- `interconnect(kind_a, kind_b, effort, flow) -> Result<Interconnection,
  CoupleError>` — a Dirac structure (shared effort, opposite flow) whose
  `interface_power` is `0` exactly (power-conserving by construction); refuses
  incompatible ports. `interface_power(&[Port])` = `Σ effort·flow`.
- `EnergyAudit` — `record`, `max_generation`, `is_passive(tol)`: the legacy
  `is_passive` name checks only caller-supplied scalar interface imbalance at
  each recorded exchange. A nonzero balance is a bug alarm, not a proof of
  whole-system or closed-window passivity.
- `AitkenRelaxation::new(omega_init, omega_max)` + `next_omega(residual)` — the
  scalar Δ² dynamic relaxation factor, magnitude-capped.
- `iterate_fixed_relaxation` / `iterate_aitken` — the added-mass interface
  fixed point under fixed vs Aitken relaxation → `FsiResult { converged, steps,
  solution, final_residual }`.

### `vibroacoustic` (bead frankensim-fsim-vibroacoustic-wgkq7)

The crate's first numerical engine: modal structure x rigid-wall-cavity
x exterior-radiation coupling, frequency domain under `e^{-i omega t}`
(matching `fs_bem::helmholtz`), hysteretic stiffness `k (1 - i eta)`.
Sign convention: positive structural deflection points AWAY from the
cavity (composes with the BEM's outward panel velocities).

- `StructuralModes` (mass-normalized, `fs_modal` convention) /
  `CavityModes` / `AcousticMedium` — the data carriers; any cavity
  basis plugs in, `rectangular_cavity_modes` (analytic, with a
  sufficiency-PROVEN enumeration box — an isotropic heuristic silently
  dropped axial modes of elongated cavities, executed regression) and
  `helmholtz_resonator_mode` (lumped lowest mode, `2 (8/3 pi) a` end
  correction; the psi = 1, Lambda = V presentation is EXACT, proven by
  the 3-state reduction) ship.
- `AcousticMedium::try_from_material` consumes the shared `FluidAcoustics`
  view, admitting original gas/mixture, source-card liquid and IF97 states
  without an ideal-gas replacement. The same density/sound-speed snapshot
  supplies cavity frequencies, reciprocal bulk loading and the exterior BEM
  medium. The casebook checks sealed-piston frequencies, displacement, pressure
  and power against the independent bulk-spring law for air, liquid water and
  dense water. Material phase and source validity remain with the original
  owner; liquid attenuation, mean flow and cavitation are not inferred.
- `assemble_coupling` — area-integral `C_rq = INT phi_r psi_q dA`,
  pinned against closed-form sin x cos overlap integrals.
- `project_radiation_impedance` — panel-space `Z` (`p = Z v`) reduced
  to modal `Zm` with area weighting on the force side only.
- `VibroacousticModel::with_structural_viscosity` — attaches the existing
  `ModalPlateViscosity` PSD owner in the same mass-normalized modal basis.
  The full signed matrix enters the coupled structural dynamic stiffness as
  `-i omega C`; passive-channel squares supply `1/2 Re(v^H C v)` to the
  structural power account. Cross-mode damping and undamped nullspaces are
  retained, including each principal block used by the truncation solve.
  Scalar structural hysteretic loss must be zero when this matrix owner is
  attached; cavity loss and exterior radiation retain separate accounts.
  No basis reprojection or frequency-dependent relaxation law is inferred.
  `src/vibroacoustic/viscosity_tests.rs` checks an independent analytic
  two-mode/cavity response, integrated cycle work with radiation, a relative
  dashpot's common-motion nullspace, signed-basis invariance, and truncation.
- `VibroacousticModel::with_cavity_relaxation` — binds the original gas's
  molecular reservoirs to the cavity through `fs-material::relaxation`'s
  complex bulk modulus `K(omega)`. Density and equilibrium sound speed must
  match the cavity basis exactly; scalar cavity loss must be zero. The pressure
  mass term becomes `omega^2 rho0 c0^2/K(omega)`, with the original reciprocal
  inertial boundary forcing retained. Molecular heat is independently evaluated
  as `omega Lambda Im(1/K) |p|^2 / 2`, including uniform sealed-cavity compression.
  This is a homogeneous, quiescent molecular constitutive law. Shear viscosity,
  thermal diffusion, wall boundary layers and mean flow are not inferred from a
  propagation attenuation coefficient. Per-frequency caloric work is bounded by
  the original reservoir count and performed once outside the modal loops.
  `src/vibroacoustic/gas_tests.rs` checks a sealed gas spring, independent Schur
  elimination, integrated pressure/strain cycle work, reciprocity, modal
  rescaling/truncation, explicit humid-air species and stale-gas/double-loss
  refusal. The undamped eigenfrequency diagnostic remains the equilibrium-gas
  reference; loaded dispersive resonances require the FRF.

- `VibroacousticModel::frf` / `frf_truncated` — dense complex coupled
  solve; every response carries the complete power breakdown with BOTH
  audit residuals (`input = structural + cavity + radiated`,
  `interface + acoustic_input = cavity`). The force-only paths retain zero
  acoustic input and use the same original coupled operator.
- `frf_with_sources` / `frf_truncated_with_sources` admit original-basis
  volume-flow projections `Q_q = integral(psi_q s dV) + integral(psi_q v_in dA)`
  with inward injection positive. The pressure row receives
  `-i omega rho0 c0^2 Q_q/Lambda_q`; the same source supplies signed work
  `Re(a_q conj(Q_q))/2`. Pressure backreaction, full modal structural viscosity,
  molecular compressibility and exterior radiation remain in the same solve.
  Source vectors retain full-basis ordering under truncation. A driven structural
  wall must not also appear as a prescribed cavity source. This admits imposed
  linear acoustic excitation; it does not infer jet physics or source impedance.
  `project_cavity_boundary_flow` maps an explicit inward velocity through the
  original cavity surface quadrature, retaining signed cancellations and the
  equality between integrated surface work and modal source work.
  The source tests compare sealed-cavity compliance with an independent closed
  form, integrate pressure/flow work over a cycle, and check reciprocal mixed
  transfers, coherent superposition, signed normalization and truncation.
- `frf_with_convergence` — truncation error OBSERVED (full vs halved
  bases) and refused (`FS-COUPLE-VIBRO-TRUNCATION-NOT-CONVERGED`)
  above the caller tolerance; near-resonance rows are genuinely
  truncation-sensitive (measured 0.28 on the casebook fixture) and the
  casebook exposes rather than hides that.
  `frf_with_sources_and_convergence` applies the same check to the total work
  of the prescribed force and volume sources. This scalar comparison does not
  certify pointwise pressure convergence or a quiet transfer-function notch.
- `undamped_natural_frequencies` — exact `x = omega^2` linearization
  (unit-triangular `B`), pinned by the two-oscillator closed form at
  ~1e-16 across a five-decade coupling sweep, a first-order added-mass
  falsifier (quartering 4.00), and a multi-mode independent-determinant
  oracle.
- `VibroError` — stable `FS-COUPLE-VIBRO-*` refusals.

### `modal_acoustic_time`

A re-export of `fs_modal::time` (IRP M1.1-W1a, owner decision OD10): the exact zero-order-hold modal realization, its types, caps and
semantics are specified in `crates/fs-modal/CONTRACT.md` (`### time`), and every former path `fs_couple::modal_acoustic_time::*`
still resolves. fs-couple keeps `modal_acoustic_time::broadband` (the opt-in broadband receiver observation; it needs `fs-vfit`,
`fs-phs` and `broadband_radiation`, so it cannot live in L1). Consumer facts that stay here: `AcousticRealizeError::Modal` wraps
`ModalAcousticTimeError`; the bowed runtime's set-valued contact uses the free-motion plus unit-force endpoint response, cached per
bow station (see the Bow bullet under the string paths); the plate and bowed-string paths call the single-mode steps
`advance_exact_zoh` and `ModalAcousticTimeModel::advance_exact_zoh_step`, which are public in fs-modal for that reason.

### `render`

Block render API (music bead `frankensim-music-v8-root-3ez8g.2.1`): the
callback-shaped hosting layer for performance images. `RenderContext`
owns admitted voices (`ReedBoreVoice` — the `realize_reed_bore` physics
restructured for blockwise stepping, with the one-shot realizer now a
thin wrapper over it so the two paths share one loop body;
`ModalStringVoice` — the exact-ZOH modal runtime hosted verbatim) and
advances them one block at a time into a pre-sized scratch buffer.
INVARIANTS: a block boundary performs no arithmetic, so block size is a
pure loop-partition choice — the one-shot path and every block partition
are bitwise-identical (tested at 64/571/full); typed `ControlDelta`s
apply only BETWEEN blocks, transactionally.
`SetModalForce` changes one mass-normalized held force [N/sqrt(kg)] in an
admitted modal voice, preserving displacement, velocity and its physical clock.
Zero force releases the input without silencing or resetting the resonator.
All voice/mode indices and force values in a batch are checked before any
input or control log changes. This exposes generalized excitation, not a
spatial pluck/contact law or an automatic score scheduler. Controls are logged as
`ControlRecord`s carrying their D17 lift description (empty = pure input
move); cancellation (`render_under_gate` over `fs_exec::CancelGate`) is
polled only at block boundaries with drain semantics — an in-flight
block completes, the rendered prefix is a whole number of blocks, and a
resumed context continues bitwise-identically (tested). ALLOCATION:
the massless-reed voice with an empty plate bank is allocation-free per
block (counting-allocator gate); the massive-reed lay path
(`dissipative_modal_forces` returns a `Vec` per sample) and the modal
voice (`step` builds a per-sample energy frame) are DISCLOSED allocating
voices — fusion candidates (bead 3ez8g.15), not silently admitted.
Determinism: one-host bitwise, inherited from the hosted kernels.
No-claims: no device/OS audio host (Franken-only policy; live output
would be a quarantined adapter decision); no image hopping yet (D17
lifts beyond input parameters land with the articulation/track beads);
refusal mid-block permanently poisons the context (mid-sample state is not
rewound). The failed block's output is partial and must be discarded.
Subsequent block or control calls return `RenderError::Poisoned` before
touching output, voices, or control logs; reconstruction is required.
Request sizing and control-admission refusals do not poison a healthy context.
`ScheduledRenderer` admits immutable `ScheduledControl` schedules.
`ScheduledRenderer::from_pressure_gestures` connects complete pressure tracks
directly to hosted voices, requiring their actual sample periods to match the
declared audio rate. It validates destination bindings even for an empty horizon,
requires context sample zero, and preserves supplied physical state. Compilation
and all admission finish before controls or simulation state advance.
The renderer
splits caller blocks exactly before their integer sample indices. Equal-index
events retain input order; end-boundary events wait for the next nonempty block.
All controls are validated and log capacity reserved before stepping. The
constructor refuses invalid schedules or budgets; no resampling or
score-to-excitation compilation is claimed. Its cancellation gate polls before each
caller block and consumes no events on cancellation. G3 tests compare scheduled
force changes, release and re-excitation against direct modal stepping across
multiple block partitions, with cancellation/resume at callback boundaries.
`pressure_gesture_controls` binds a selected blowing-pressure gesture track to
that scheduler. It samples existing gesture ramps at the declared control rate,
holds each value, and places ticks at the first audio sample at or after their
time using integer clock conversion. Admission bounds evaluated ticks and storage;
the caller supplies the destination voice's audio rate and a half-open horizon.
G3 checks a nonintegral clock ratio and pressure ramp against per-sample control
application across callback partitions. Other gesture targets and continuous-time
interpolation between control ticks are not implemented by this binding.
Boundary facts a consumer must
know, both discovered by this bead's battery: the characteristic-line
realization refuses `UnflangedOpen`/`FlangedOpen` terminations whose
Nyquist `ka` exceeds 1 (at 48 kHz that is every bore wider than
~2.27 mm radius — the measured-load `MouthLoad::Tabulated` lane
is the lift), and the Zwikker–Kosten Bessel path's former
`r_v ≈ 160` continued-fraction ceiling was an iteration-budget bug fixed
in `fs-phs` (budget now scales with `|z|`).

The scalar characteristic-line constructors, `realize_line`, `impedance_line`,
regional/material-wall lines and all duct/field reed-voice constructors accept
the shared `MouthLoad` selection. Existing `Termination` arguments retain the
analytic path. A tabulated load requires an explicit `TabulatedMouth`: the
mouth plane/radius and exterior gas must match, and its real DC impedance is
supplied separately rather than extrapolated. Every positive-frequency bin
comes from the table, with no analytic splice or out-of-band fallback. The
complete declared FFT/ladder band is admitted before realization. Raw terminal
data alone cannot supply the bore's static behavior: the gas/wall propagation
must also admit the existing 1e-3 rad/s quasi-static probe. A material wall
whose source band excludes that probe explicitly refuses, even when its
positive FFT band is covered; no static wall response is extrapolated. Raw terminal
passivity is checked before propagation, so bore loss cannot conceal an active
mouth; characteristic FIRs also retain the 4x-oversampled passivity refusal.
The existing Bessel propagation, FIR/ladder and reed stepper remain the numerical
owners. Scalar lines refuse `Matched`, which is a multimodal continuation.
The full-band regressions compare characteristic reflection and mouth-pressure
transfer with frequency-domain TMM at 48/96 kHz, plus uniform-regional identity,
changed-load sensitivity, physical-reference/band refusals and active-load
mutations. They use an explicit resistive control, not a validated exterior
radiation bake; table convergence, complete-instrument realism and multimodal
runtime adoption are separate work.

### Aperture-junction solver modes (bead frankensim-2s4i5)

IRP M5.5 adds `render::ladder::RealizedResolution` and
`ReedBoreVoice::{new_with_fft_length, realized_resolution,
set_junction_step_tolerance, junction_residual_max, initial_ir_tail_energy}`.
The constructor shares the existing `with_line` admission/initialization; it
only passes the explicit request to `characteristic_line_dense`. Actual FIR
history length exposes the 8192 clamp. Strict mode has **no** pressure-step
tolerance: it still terminates on an exact zero or adjacent binary64 endpoints.
The massless FastNewton seam admits an explicit relative tolerance in (0,1],
keeps 1e-9 as the existing default, and preserves every corner/fallback rule.
Massive reeds refuse this tolerance seam rather than reporting a tolerance the
coupled structural solver did not use. The observed volume-flow defect is
recorded in m3/s, and the initialized periodic FIR's last-quarter tap norm is a
dimensionless diagnostic. Neither bounds pitch, onset, or partial-level error.

`tests/convergence_ladder.rs` binds the existing M-CL0 fixture to fs-ladder,
with at least three measured sample-rate, FFT, and Newton-tolerance rungs and
the joint group at fixed physical IR duration. The legacy observer's backward
jet-force difference supplies the declared smooth-signal first-order hypothesis;
this is not a global-order assertion for the oscillator. The fixture measures
raw pitch/partials and a **finite-horizon onset** (0.25s, 1e-3 Pa RMS, nondecaying
last-quarter RMS, 8 bisections in the declared 0..4800 Pa bracket). An unbracketed
onset refuses. This is not an asymptotic bifurcation/extinction measurement.
G1/G3/G5 report construction is exercised, but Q2 accuracy authority remains
unavailable where no conditioning map supplies a bound in the acoustic QoI's
units. Spreads and signed finest-reference joint interactions remain data.
`scripts/e2e/music_ladder.sh` runs a prebuilt native test executable and retains
JSONL, a raw-scaled PCM render/provenance sidecar and its fs-blake3-bound report
in a fresh directory. It performs no build and never overwrites evidence.

Both modes share the closed-slit characteristic constraint
`(1-r0) p_plus - p_minus_history = Zc Ubody`. The reflection denominator
is not clamped, and closure is evaluated at the body-loaded bore pressure;
body motion can reopen the slit. Singular or nonfinite candidates do not
qualify for this shortcut. G1 junction tests cover signed body flow,
reflected incoming waves, and reopening.

The reed aperture junction has two solver modes on
`ReedBoreVoice::set_solver_mode`. `Strict` (default) is the
deterministic path: a local sign-changing bracket, up to 32 symmetric
bracket expansions if needed, then up to 2100 bisections (covering the
binary64 exponent range, including subnormals). The former 96-step bound
could refuse near-zero roots solely because their pressure scale was far
below the initial bracket width. It stops at an
exact residual zero or adjacent floating-point pressure endpoints,
returning the endpoint with smaller flow residual. Nonfinite residuals,
failure to bracket, or exhausted refinement return errors; a grid guess
is never returned as a solution. This removes the former dimensionally
invalid pressure-scaled flow tolerance. Root selection is deterministic,
but multiple quasistatic roots and the Bernoulli dead zone remain model
limitations; this is not a uniqueness or exact-arithmetic certificate.
G1 tests independently reconstruct volume-flow balance across body flow,
reflection, and large incoming pressure, and test refusal without a root.
`FastNewton` is a DECLARED FAST MODE (the
fs-rand ziggurat precedent): DAMPED island Newton on the analytic
Jacobian of the shared Bernoulli-plus-wave-junction residual — each
step is halved up to four times until the residual magnitude strictly
decreases, and a sample that never improves, hits a sqrt kink or
non-finite state, exhausts its iteration budget, or shows a vanishing
slope is handed untouched to the strict path and counted in
`FastSolveStats`. Convergence requires an UNDOUBLED full step inside
`NEWTON_STEP_TOL`, so a still-descending line search never reads as
a converged root. It is NOT bitwise-equal to the strict path and must
never become the default until it earns the same proof. Receipts (bead
2s4i5): the junction-level battery lives inline (`fast_mode_tests`) —
a 64-point open/interior/closing stimulus sweep checks both modes against
a volume-flow scale `w H sqrt(2 Pc/rho)` and checks pressure agreement
within `1e-8 * (1 Pa + |p_strict|)`. The integration lane
(`tests/reed_newton_fusion.rs`) owns the fallback-hit-rate receipt on
a full nominal render (bounded structurally at 0.5; measured 12-22%
warm-seeded) and before/after budget rows stamped with build profile —
debug rows are diagnostics forever, per the budget-lane doctrine.
Recorded finding: end-to-end rendered waveforms are NOT comparable
between modes beyond onset — the self-oscillating loop amplifies
per-sample root differences into phase-divergent (O(kPa)) waveforms,
so waveform equality is never a claim either mode makes.

### Material-driven viscous slit flow

`bernoulli_aperture::viscous::ViscousSlit` owns the shared reduced channel law
`dp = 12 mu L U/(w h^3) + rho U|U|/(2 w^2 h^2)`. `from_gas` retains one
coherent `GasState`, derives viscosity from it, and refuses a different density
at evaluation. `PairedAperture` and `ModalAperture` can select that law from
their retained material gas plus the physical channel length. `DynamicAperture`
accepts the same source-bound law; a plate-valve performance's optional
`slit_length` record uses its actual regional inlet gas, not the exterior gas.

All three mechanical owners use this flow inside their same-sample pressure,
contact and swept-volume solve and account for the channel loss as `dp * U`.
Independent blades and signed modes retain their own conjugate pressure areas;
the shared slit width is counted once. Spatial single/symmetric/independent
plate closures integrate the same law analytically along each positive P1 edge
segment, with an analytic pressure tangent. Closed portions contribute zero.
An area-equivalent height is never substituted into the cubic viscous term.
Spatial pressure condensation uses safeguarded Newton with a passive bracket,
96 evaluations and a 32-epsilon relative pressure-balance target; failure
refuses the trial before any mechanical, material-memory or wave state advances.
Law attachment is allowed only before the first accepted step.

G1/G3 tests in `viscous_aperture`, `modal_aperture::viscosity`, the paired spatial
module and the constitutive owner cover flow/work closure, material-gas swaps,
the cubic-gap laminar limit, independent point quadrature, partial closure,
edge subdivision, paired/modal equivalence and full modal damping. This remains
a fixed-state, parallel-channel, quasi-steady model with one common streamwise
length. It does not resolve transverse/entrance flow, channel inertance,
pressure recovery, gas evolution, compressible flow or measured instrument
accuracy; it is not the plan's resolved-flow acceptance rung.

### `music_render` (binary)

`ReedPerformance` accepts optional `slit_length LENGTH_M` immediately after its
reed record. `ReedBoreVoice::set_viscous_slit_length` binds the shared
`bernoulli_aperture::viscous::ViscousSlit` to the exact source inlet `GasState`,
including regional mixtures and axial fields. Both quasistatic and massive
reed/contact branches solve this constitutive law in the original same-sample
characteristic junction; no output filter, alternate contact owner or new
integrator is introduced. Omission retains the inviscid arithmetic. Selection
after rendering refuses. The massless viscous branch uses strict refinement;
fast mode counts this fallback and exposes no Newton stopping-tolerance claim.
G1/G3/G5 regressions in `render::scattering`, `reed_bore::fast_mode_tests` and
`render::schedule::reed::viscous_tests` cover the signed pressure/flow equation,
Poiseuille/closure limits, ten beating/reopening cycles with a pressure-work
account, inlet-material binding, and callback invariance. Quasisteady planar
flow does not establish unsteady channel, turbulence, or instrument realism.

Staple pressure recovery (`staple_recovery C_P S_C_M2`, MB 5047, 2026-11-05). A double reed's channel
widens from the constriction of area `S_c` into the conical staple, and a fraction `C_P` of the dynamic
head at the constriction comes back as pressure (Almeida, Vergez and Caussé 2007, J. Acoust. Soc. Am.
121(1) 536-546, arXiv physics/0607011v1 pp. 8-9, eqs. 14-17). An optional reed-card record
`staple_recovery C_P S_C_M2` right after `slit_length` (so it needs that row, a massive reed and a fresh
voice) selects it through `ReedBoreVoice::set_staple_recovery`. The flow through the opening follows
the viscous slit law on the constriction drop `(dp)_c = p_m - p_c`; the blade is driven by `(dp)_c`
(eq. 16), not by the total mouth-to-bore drop (`JunctionPressure::blade_drop_pa` in `reed_bore::contact`,
taken by the momentum balance of both propagators); and the bore input is
`p_r = p_m - (dp)_c + C_P (1/2) rho (q/S_c)^2` (eq. 17), with `q` the flow through the opening.
`ViscousSlit::characteristic_junction_with_recovery` solves the same-sample characteristic junction: for
forward flow the slit's quadratic term takes the density `rho (1 - C_P (S/S_c)^2)`. Recovery acts on
forward flow only (in reverse the staple is upstream of the jet, so there is nothing to recover) and on
the JET flow: the moving blade's swept flow and any
body flow are not carried through the staple, which is where the quasi-static measurement stops being
one. The law is singular where `C_P (S/S_c)^2 >= 1`, where it would cancel the whole jet loss; there it
refuses (`STAPLE_RECOVERY_DOMAIN`), and `set_staple_recovery` refuses a reed whose rest opening already
is there. `C_P = 0` is the original law and installs nothing: the row absent, or with `C_P = 0`, is the
voice as it was, bit for bit. The only printed numbers are one measured reed's, `C_P = 0.8` (fitted on both
pressure branches, Fig. 15) and `S_c = 4e-6 m^2` (Fig. 14): a ready-to-play natural-cane oboe reed in
humid air, masses added, no lips. A card that uses them borrows them for another reed, and says so.
Energy: inside the domain the slit and staple element dissipates `(dp)_total q = viscous + (1/2) rho v^2 q
(1 - C_P (S/S_c)^2)`, never negative. Because the blade feels `p_c` while its swept flow enters the bore
at `p_r`, the step's energy account carries one exchange term, `-(C_P (1/2) rho (q/S_c)^2) * swept flow`,
which has either sign. G1 tests: the flow maximum of the printed characteristic sits where the closed
form of eqs. 15-17 puts it (the total drop there is `x (1 - C_P S^2/S_c^2)` at `x = p_M/3`,
`S = 2 S_0/3`, computed in the test); the domain refusal by name; the junction's power balance and
passivity; the blade momentum under `(dp)_c` and the step's energy account; bit-identical absence;
callback invariance; the card rows' refusals.

The music lane's render CLI (bead ib15w): `music_render <reed|string>
<out.wav> [--seconds S] [--block N] [--full-scale-pa P]`. Pinned fixture
compositions render through the block API and encode through the ONE
pascals→PCM owner, `pcm_wav::encode_pcm16_wav` — the recorded seam
decision (beads ib15w + h7xu5.7.8): the cinematic stereo/receipt-hashed
encoder stays cinematic, no third RIFF writer exists in the music lane.
Reed runs accept `--temperature-k K` and `--ambient-pressure-pa PA`; these
feed the shared `GasState` into the existing duct/reed calculation.
Optional `--relative-humidity RH` supplies a dimensionless fraction in `[0, 1]`
to the shared `GasState::try_new_moist_air` mixture owner. Default zero follows
the exact dry-air path. Nonzero humidity changes the gas properties consumed
by the duct and reed calculation; it is not an output pitch/gain adjustment.
For nonzero humidity the owner refuses temperatures outside its saturation-fit
window and excessive vapor fractions. This does not model condensation or
moisture uptake in solids. The ambient values, relative humidity and gas
specification are recorded in the sidecar. Other fixtures refuse these
overrides. This is a stationary gas-state change, not
solid thermal expansion, heat transfer or phase change. Reed runs also accept
`--duct-length-m M` and `--duct-radius-m M` for the cylindrical duct, feeding
the declared geometry to the existing acoustic solver and recording it in the
sidecar. These dimensions must be finite and positive; other fixtures refuse
them. Optional `--duct-outlet-radius-m M` selects a linear cone when different
from `--duct-radius-m` (the inlet); omission or equal radii retains the cylinder.
Both cone endpoint radii are recorded for replay. Other fixtures refuse
them. Existing acoustic admission also applies: a positive radius can still
exceed the selected radiation model's `ka` limit and must refuse.
Geometry remains fixed during the performance. No automatic wall-material
response is implied. Gas-model validity
limits do not certify the fixed solid specimen at extreme temperatures.

The file-driven path also accepts `frankensim-flue-performance-v1`:
`music_render wind crates/fs-couple/examples/flue.performance OUT.wav --block 512`.
The fixed-order SI records explicitly supply mechanical/output clocks, the output
sample window and pascal scale/ceiling, bore and exterior gas, cylinder dimensions
and prescribed terminal reflection, mouth dimensions, jet card/gain/onset seed, an
optional opt-in jet response in place of the gain chain (`jet_response forced-bickley
SPACING WINDOW [exit kutta|measured]`, `jet_transfer nolle-price-2015 CHANNEL CHAMFER_DEPTH_M` or
`jet_transfer auvray-2014 BETA GAMMA SOURCE_TAG`, exclusive; the last is Auvray et al.'s
jet drive (JASA 136(1) 389, 2014), an in-phase receptivity `h exp(BETA W/h) v / U` delayed
by `W / (GAMMA U)`, a Bickley cut `b = 2h/5` and Verge's source length `(4/pi) sqrt(2hW)`,
on `U = sqrt(2P/rho)`, so it refuses every windway row; BETA and GAMMA must be positive,
GAMMA at most 1, and no field has a default; that paper's Table I gives BETA 0.3 and
GAMMA 0.4, validated on one recorder at W/h 4; `exit measured` is admitted only where
Re_h = U0 h / nu lies within a factor 2 of 200, the one Re_h at which de la Cuadra measured
it (Fig. 2.28); elsewhere the voice refuses),
jet limits, wave memory, fractional phase band/budget, optional viscothermal cells
and band, receiver position/resolution/band, and control compilation limits. The
suffix is the existing canonical, typed blowing-pressure gesture schedule.
`jet_card tonal-interim` explicitly selects the retained stage-I source card;
`jet_card json {CANONICAL_JSON}` supplies a complete one-line card. No filename or
instrument name selects physical dimensions, material properties or calibration.
`gas` and `exterior` accept the same legacy, named-mixture and source-backed NASA-9
records as other wind sources; `exterior same` shares the supplied bore gas. A
separate exterior must have the same static pressure. Gas changes enter jet speed,
Reynolds number, propagation, selected viscothermal loss and physical radiation.
The optional `windway none`, `windway slit LENGTH_M` or
`windway unsteady LENGTH_M TRANSVERSE_CELLS CONVECTION_CELLS` record after
`flue` selects the shared ideal, quasi-steady or unsteady source-gas windway
before performance setup. Actual channel length, model, viscosity and resolved
cell counts remain in the output sidecar. Omitted/`none` retains the ideal jet.
A finite channel is explicit physical geometry; no lip-channel length is
invented for a transverse flute.

`windway slit-developing LENGTH_M` (opt-in, pexfp.135.4) holds the jet by
`bernoulli_aperture::developing`: a straight parallel-plate slit fed through a
uniform inlet, `P = rho U^2/2 [1 + 4 f_app L/D_h]` with Shah's (1978) apparent
friction factor for developing laminar flow, as stated openly by McCreery et al.
(NURETH-15 2013, eq. 2). No constant is tuned; at large `x+` it is the
fully developed slit law plus `0.674 rho U^2/2`. The held note's pressure (the
gesture's highest blowing pressure) must lie inside the correlation's validity
or the input is refused: `x+ >= 4e-4`; Chen's low-Re inlet term
`0.64 + 38/Re` may move the held drop by at most the stated 2.4 % (Chen's term
and Shah and London's agreement statement behind the `x+` floor are secondary,
search snippet only; primary not read); an
undeveloped exit keeps `Re_delta* <= 520` (core speed and displacement thickness;
the Blasius layer's linear instability, Manneville 2016); and a developed exit
keeps `Re_Dh <= 1760` (`Re_cl = 660` of the laminar centreline speed and half gap,
above which plane Poiseuille flow sustains turbulence, Kashyap et al. 2020). The
sidecar records `x+`, `Re_Dh`, `Re_L`, the core
speed, displacement thickness and these measures at the held pressure. The
attack and release pass through the same quasi-steady law and are not covered.
The forced Bickley jet and the Nolle/Price transfer, attached after it, take its
exit at the held pressure: the law's core speed over its bulk speed, exactly the
Poiseuille 1.5 once the layers have met (see OPT-IN FORCED JET; pexfp.135.5). A
windway the law refuses never reaches them, and a transfer attached before this
windway is refused. Price's jet gain 0.145 was fitted on a Poiseuille jet and is
applied to this exit's centre speed unchanged; the transit and gain follow the
held ratio, so the attack and release are not covered. Behind either jet the
island's card speed (the onset seed, the Mach admission and the reported Reynolds
number) is that jet's centre speed at the held ratio, while the flux keeps the bulk
speed; the low-Mach domain therefore binds the centre speed, and a uniform exit's
ratio is exactly 1. Without a jet transfer the card reads only the bulk speed
unless `jet_speed card-centreline` is named. A
converging windway declared as its
exit-height slit most likely overstates the loss.

`jet_speed card-centreline` (opt-in, pexfp.135.5) after the `jet` record feeds the
jet card the exit's centre speed. By default the card reads the windway's bulk
speed, but its Strouhal number and band are defined on the centre (peak) speed of
a near-top-hat jet (`MeanJetProfile::u_centerline`). The rule is:

| windway | card speed |
| --- | --- |
| none | the bulk speed `sqrt(2P/rho)`; nothing changes |
| `slit-developing` | the inviscid core `sqrt(2P/rho)` while the layers have not met, then the Poiseuille centre `1.5 U` |
| `slit` | the Poiseuille centre `1.5 U`, refused where it exceeds `sqrt(2P/rho)` at the held pressure |

The card speed sets the convection `c_p = (c_p/U) U`, the growth band with its
band-lag tap, the onset seed and the card Reynolds number. The bulk speed alone
carries the flux, so the held flow is unchanged. The unsteady channel and the jet
transfers, which own their exit speed, refuse it. The sidecar records the law and
the held bulk and card speeds. Omission keeps every default result bit-identical.
No card selects it: under it the island's loop-phase window moves every flue card
off its pitch (pexfp.135.5 A/B), and that window is the open defect.

`jet_noise SEED GAIN ST_LO ST_HI RE_ON` (opt-in) after the `jet_limits` record adds a
measured turbulence source to the mouth. Every field is required and none is
defaulted: the card names the published source it uses. `U` is Verge's jet
velocity `U = sqrt(2 P_foot / rho)` (Verge 1995, eq. 2.12 definition) of the
instantaneous blowing pressure `P_foot` upstream of the windway, for all three uses
below; it is not the windway's bulk or centre speed that the loop's jet has, so at a
fixed `P_foot` a card with another windway loss draws the same noise (a test pins
this). The amplitude is `GAIN * rho * U^2 / 2 = GAIN * P_foot` in rms INSIDE the
band, the gate is `Re_h = U h / nu >= RE_ON`, and the band is the Strouhal band
`St_h = f h / U` between `ST_LO` and `ST_HI`, recomputed every sample from the live
`U`, so the spectrum slides with the pressure and the rms follows `U^2`. The stream is
a unit Gaussian through one 2-pole Butterworth band-pass (one biquad, Verge's own
recipe, sec. 6.5 p. 196, which he never ran, p. 198) whose -3 dB edges are the band's
ends, prewarped so they are exact at the mechanics clock, and scaled by the exact
variance the section has inside the band (16-point Gauss-Legendre), so `GAIN` is the
in-band rms; outside the band the skirts fall at 6 dB/octave, like the pedestal of
Verge's Fig. 5.15a (p. 160), and the whole stream's rms is about 1.4 times the
in-band rms. A band that reaches the mechanics clock's Nyquist frequency keeps only
its part below it, at the whole band's density, and a band wholly above it is silent
(the rule the voice's and the lip's rows need; no flue band reaches Nyquist, and
below it every sample is bit for bit the band-pass alone). The stream is keyed by
`SEED` and its own kernel id
(`FLUE_JET_NOISE_KERNEL`), not the island's onset-seed stream: exactly one normal is
drawn per sample whether the gate is open or shut, so a gate that never opens leaves
the loop bit-identical to the card without the row. The gate returns zero (and
clears the filter state) while `Re_h < RE_ON` or the jet is not flowing. The
pressure enters in series in the mouth solve, `radiation_history - p`, with the same
value in all four trial solves of a sample; its work `p q dt` is booked in the
source-work ledger and the ledger still closes. Omission keeps every result bit-identical. The sidecar records
the five fields. No-claim: a single broadband source in series at the mouth; the
gain and band belong to the cited measurement (Verge 1995 sec. 5.5 for a recorder-like
pipe), are valid only near the Reynolds number and W/h where it was measured, and
nothing here shapes the spectrum beyond one section's skirts or tunes a gain to a
recording. Until 2026-11-03 the band was a difference of two one-pole low-passes
normalised on the whole band, which put 0.26 of its power inside the band and so
delivered `GAIN` 5.8 dB low there; that filter (`BandedGaussian`) remains, unchanged
and used by no row, as the tests' planted reference. The voice's and the brass lip's
jet rows draw through this row's band-pass, under their own kernel ids.

`bernoulli_aperture::slit_channel` resolves fixed, fully developed planar
Newtonian momentum from the retained `GasState` density and viscosity. Its
2..=4096 cell-centred finite volumes share a no-slip wall diffusion operator;
backward Euler evolves the transverse profile from rest under signed pressure.
An optional uniform kinetic-head exit is solved with the mean flow in the same
implicit update. Kinetic storage, pressure work, physical viscous heat, exit
work and numerical dissipation remain separate. Standalone refusal leaves the
accepted channel state intact. No fitted relaxation poles or alternate gas
formulas enter this owner.

`FlueVoice::with_unsteady_windway` connects that accepted flow to jet speed,
source volume, Reynolds coverage and the original card's growth law. Pressure
controls accelerate the fluid and release retains its momentum. A separate
2..=4096 cell semi-Lagrangian convection field transports the existing feedback
across the actual mouth cut-up at the evolving speed, including stop and
startup, without a minimum-speed history clamp. The labium feedback is read
at the frozen path's own delay (pexfp.320, 1718d79fb): the transit less the
band filters' low-frequency lag `2/omega_c + 1/omega_h = 16 theta / U`, with
its one-sample floor and sub-sample rule, as the distance `c_p tau` along the
field, so `16 theta c_p / U` upstream of the cut-up, independent of `U`. The
read interpolates linearly between the bracketing nodes (between the current
mouth velocity and the first node below it), the transport step's own
interpolant. A sub-sample transit and a jet at rest read at the labium. Before
this the transported loop was longer by `16 theta / U`, which is 88 to 114
degrees of jet phase at f0 on recorder a1/f1 and whistle D5. The voice's
original failure poisoning applies to the combined sample. Channel and acoustic mouth work
ledgers retain their distinct ports; numerical damping is not material heat.
Characteristics entering during a sample interpolate the accepted old boundary
velocity and the current solved mouth velocity at their physical entry time.
The affine current coefficient participates in the existing joint jet/mouth
equation; interior characteristics use the preceding spatial field. Zero speed
preserves interior parcels and retains the accepted boundary trace for restart.
Channel momentum, convection and bulk speed remain private candidates until
the mouth and wave step accepts. A rejected trial preserves their preceding
accepted state under the existing poisoned-voice policy. Convection uses a
piecewise-constant accepted bulk speed during each source sample. Focused
characteristic, source-history and gas-to-PCM regressions are authored; syntax
parsing and an independent characteristic/refinement oracle passed. Native
execution remains required.
Single-tube and regional-gas branched file sources use this same builder and
the same physical receiver. Independent steady Poiseuille and transient
sine-series regressions exercise space/time refinement, signed release,
reversal, gas interventions and callback partitioning; native execution is
required before treating them as passing evidence.

These are explicit numerical axes for a fixed, incompressible, fully developed
planar channel and the original estimated jet card. Side-wall/entrance/axial
profile development, compressibility, thermal transport, normal wall motion and
turbulence remain outside the channel. The exit uses uniform bulk kinetic
head despite the resolved internal profile. This does not establish a resolved
jet, moving-reed channel or experimentally validated instrument claim.

This loader constructs the existing `FlueVoice::from_fractional_network` and
`FlueOutletPerformance::with_multirate_radiation`. Every mechanical outlet-flow
sample reaches the shared differentiating Rayleigh receiver before the sole PCM
owner. The file declares 48 kHz CLI output and an integer mechanical/output ratio;
`--decimate` refuses a second observer. The sidecar retains both clocks, original
input/card identities, source gas coefficients and the actual receiver delay.
This is an **EstimateOnly** composition of a tonal jet and stationary rigid
cylindrical pipe with a prescribed reflection and one-way infinite-baffle outlet
radiation. The selected card keeps its narrower laboratory scope. It supplies no
resolved turbulence, matched radiation feedback, measured flute geometry, wall
material response or calibrated SPL claim; `--acceptance` refuses this source.

### `bakeoff`

Bake-off receipt schema and protocol for the music-program claims registry
(bead `frankensim-music-v8-root-3ez8g.1.4`). `BakeoffReceipt` records two
contender images run on one shared-card fixture against a caller-supplied
reference QoI map, with logical budgets (states, steps, solver iterations —
never wall-clock, so bytes are host-stable), observed failure modes, and a
reviewed outcome. `BakeoffOutcome` has exactly three variants —
keep-both, keep-for-subset, refuse-newcomer — and structurally cannot
express deleting an admitted image (doctrine D21: menus, not winners).
Receipts serialize to a canonical line-oriented byte encoding
(`frankensim-bakeoff-receipt-v1`, strict round-trip decoder) and are
content-addressed via `fs-blake3` under
`org.frankensim.fs-couple.bakeoff-receipt.v1`. The harness MEASURES;
the outcome is a reviewed judgment passed in by the caller — residual
tables cannot adjudicate claim scopes alone. The executed
`tests/bakeoff_string.rs` fixture (modal-ZOH vs the Gonzalez-stepped
pHS modal bank, analytic damped-oscillator reference) golden-pins its
receipt at `tests/receipts/string-modal-zoh-vs-phs-bank.bakeoff`; the
pHS contender's phase error is gated against the implicit-midpoint
dispersion model `(2/dt)·atan(omega·dt/2)` rather than a blanket state
band, because pointwise residuals near zero-crossings are unboundedly
phase-sensitive. No-claim: a receipt proves the comparison ran and what
it measured, never that a claim scope is correct — scope promotion is
the per-track gates beads' review job through the instrument-claims
registry.

### `broadband_radiation`

`broadband_radiation::modal` and `modal_acoustic_time::broadband` add an opt-in
one-way pressure receiver over complete mass-normalized modal velocity fields.
Computed training/independent held-out samples use exp(-i omega t), with explicit
producer identity, SI units, pure delay and producer discretization estimate.
Continuous vector fits are conjugated for the fs-vfit owner; proper stable
state-space filters join each original mode in an augmented held-force ZOH
exponential, reusing fs-phs's checked deterministic exponential. There is no
Tustin warp or velocity difference. The original modal owner still advances
mechanical states and work/loss records, preserving those bits exactly. Its
legacy pressure coefficients are zeroed only in the observer's private clone,
so the admitted receiver owns the pressure ceiling. A failed step publishes no
filter/delay/mechanical state. Block render failures poison their stream wrapper.

`BroadbandModalObserver::prepare_linear_forcing` additionally prepares exact
first-order-held force observation at the original sample clock. The continuous
receiver is augmented with one force-increment coordinate in normalized time;
only its extra filter moment column is retained. `step_linear(start, end)` then
advances the receiver along the same force ramp used by the shared modal owner
for motion, signed input work and physical viscous heat. The warm step reuses
allocated buffers, with no fitting, exponential preparation or transcendentals.
All pressure/mechanical gates precede publication of mechanics, filter memory
and delay together. Continuous filters retain bitwise held-step compatibility
when endpoints coincide. Constant real velocity gains admit free/soft motion;
constant complex narrow-band quadrature explicitly refuses arbitrary ramps.
This remains a one-way receiver and adds no mechanical load or damping heat.
The existing sampled-delay approximation and producer/fit error scope remain.
Six focused unit regressions cover independent polynomial convolution and
physical ODE/work/heat integration, time subdivision, rollback, held-step
compatibility and real/complex feedthrough boundaries.

The pure delay uses exact integer taps or an explicit centered Blackman/sinc
fractional kernel; all taps must be causal within the physical delay. History
starts at zero. This approximates a bandlimited sampled pressure history, not
continuous retiming of all held-force spectral images. Held-out error includes
the actual kernel and retains EstimateOnly authority; it is not a continuum
certificate. Producer error remains separate. Constant samples with zero delay
use the exact legacy equilibrium-relative quadrature operation order.

Caps: all 4096 modal inputs are retained in <=256-input cold batches, each grid
has <=512 frequencies, <=4194304 transfer cells total, per-mode order<=128,
iterations<=128, delay-kernel half-length<=4096, delay history<=1048576 samples,
and total augmented matrices<=4194304 entries. No mode/pole truncation or implicit
budget enlargement is allowed. Host cancellation is at bounded pressure-block
boundaries; no hard-real-time, continuum passivity or load-model claim follows.

`BroadbandModalObserver::{prepare_harmonic,step_harmonic}` and
`{prepare_tones,step_tones}` follow continuous cosine/sine modal forces through
the original fitted receiver, including the mechanical startup transient.
Forces are interval-local peak coefficients in exp(-i Omega t) convention;
multi-tone arrays are mode-major in the caller's original frequency order.
Prepared plans bind the actual receiver/filter identity and original sample
period. Phase advances only after a successful observation. The original
fs-modal owners still produce the mechanical states, signed input work and
independent positive viscous heat, retaining cross-frequency power terms.

For each distinct frequency, the same deterministic exponential augments
[scaled displacement, velocity, continuous filter, cosine force, sine force].
Only its two filter-forcing columns are retained; all tones share the existing
homogeneous filter transition. Thus warm steps allocate nothing, perform no
matrix exponential or transfer refit, and retain each mode's homogeneous
motion once. Exact duplicate frequencies use original-order compensated force
sums before filter application; near frequencies remain distinct and retain
beats. Cold observation preparation keeps the aggregate augmented-matrix cap;
the mechanical owner separately admits the explicit `ModalToneBudget` and
its 32-tone hard limit. The original Nyquist, fit, sampled-delay, producer-error
and one-way-load boundaries remain unchanged.

Mechanical states, the last accepted work/heat frame, filter history and delay
commit together after all gates. A spare preallocated workspace keeps the last
accepted frame intact even when mechanics refuses after acoustic staging.
Held, linear-ramp and harmonic steps share this commit boundary. Real constant
velocity gains remain valid, including free coordinates; complex constant
legacy quadrature refuses the continuous-forcing paths because it does not
define their causal broadband startup response. G0/G1/G3 checks cover the
public sampled-fit construction, independent analytical transient convolution,
true-frequency pressure and sampled delay, bit-identical original mechanics
and work/heat, DC agreement, duplicate-tone cancellation, and warm refusal /
retry across harmonic and ramp steps. These are finite-model numerical checks,
not independent experimental source or continuum-radiation qualification.


The v1 modal-file optional `broadband` block reuses its original parser and
ProjectedForceProgram; v2-v5 coupled/nonlinear trajectories refuse this
independent-mode receiver. CLI sidecars retain every producer/fit/error/delay
field, and the legacy path remains the default. The native G0/G1 battery is
`tests/broadband_modal_transfer.rs`: independent analytic held-input sphere
endpoint response at five off-resonance drives (1e-6 amplitude/radian bands),
constant compatibility, bitwise mechanical energy/work/loss, all-mode cap,
held-out refusals, integer delay, and three-rung delay/pole/grid ladders. The
pole/grid axes use the independent nonrational Bessel/Struve piston owner:
the delay-stripped sphere is already first-order rational and cannot establish
pole-count liveness. Numerical/native qualification is pending. The audible
replay/refusal sub-lane is `scripts/e2e/music_modal_broadband.sh`, requiring the
native qualification log and a supplied binary; it makes no instrument-realism
or general source-mortar/loading claim. M-BAR-Q2/M-GT-Q2/M-VNP-Q2 must supply their
own computed mode-basis samples explicitly before adopting this observer.

Generic offline bridge from solver-neutral sampled scalar-input radiation transfers to
real-tesseral SH filters. Radiation solvers supply complex-SH training rows, direct
held-out fields, source semantics, and diagnostics; `fs-couple` does not depend on them.
It converts `exp(-i omega t)` to `fs-vfit`'s `s=+i omega`, fits prewarped abscissae, and
applies unprewarped Tustin. Direct disjoint held-out fields gate artifacts; stable
proper fits may retain caller-admitted constant feedthrough but always forbid `s*e`.
The runtime owns visible per-filter state, superposes in fixed input order, and refuses
whole samples transactionally. EstimateOnly covers the fixed reference before `1/r`;
moving-source, FW-H, near-field, feedback, room, head, and passivity remain unclaimed.

### Layer-specific material memory in the reciprocal string bank

`with_material_layered_string_bending_relaxation` selects an independent
equilibrium-modulus role and ordered Maxwell-arm property pairs for every
original core/coating layer. Layers may have different names and arm counts;
an explicitly empty spectrum represents an elastic layer without invented loss
data. The existing uniform-selector API lowers through this same path.
`fs-material::section::relaxation` now owns source admission and section
integration instead of duplicating them in the acoustic coupling layer.

The bank retains its existing mass-normalized curvature projection, per-arm
history, reciprocal body/bridge loading and exact propagation. Explicit
Kelvin-Voigt viscosity remains a separate parallel dashpot; supplied Maxwell
increments must describe additional relaxing stiffness. Source frequency bands,
memory-state budgets, cold binding and sample-clock admission still apply.
Source images expose normalized `layered_selections()` plus the original
sections, layer spectra and channels. Focused G3 regressions cover bitwise
uniform/explicit layer equality and a two-arm coating/elastic-core substitution
through actual pressure and sympathetic body response, retaining the original
energy-balance assertions. These are numerical coupling checks, not physical
instrument qualification. Helical winding, axial/torsional memory, evolving
temperature and finite-strain material history remain outside this binding.

### Source-defined bridges with omitted string compliance

`ReciprocalStringBank::with_material_bridge_interfaces` and
`with_residual_flexibility` compose in either order for ordinary sourced
spring/dashpot joints. Original material identities, selected state, validity
domains, area and physical k/c remain in the source; the existing assembly
condenses the retained connection once to k/(1+kG) and c/(1+kG)^2. Later
source reconstruction retains that distinction. Binding a source to an already
corrected live system preserves its motion and pressure history. The native
`material_strings` entry validates the final condensed frequency bound.
The spring reduction is static; the dashpot is first order in c G omega.
Constant-loss anchors require a different series law and refuse in either
binding order. Focused G3 tests compare unchanged trajectories under binding
order and source attachment, then vary the joint data through actual reciprocal
pressure. These do not establish broadband continuum or instrument fidelity.

### Material string bank rotational bridges

`bowed_string::runtime::reciprocal::bank` admits explicit small-angle
`StringBodyRotationConnection` laws between a sourced torsional string mode
and a body rotation. The same polar inertia and GJ/shear-viscosity sources that
drive bow surface motion supply the string angle map. A body rotation can come
from two retained bridge footprints in the original material plate basis and
their physical separation; this is a secant rotation, exact for rigid rocking,
with a spatial approximation for a deforming plate. It is not an independently
tuned body mode. Raw angular modal maps retain their producer's authority.

The existing coupled-modal owner receives the exact work-conjugate conversion
`x = r theta`, `k_x = k_theta/r²`, `c_x = c_theta/r²`, using the sourced outer
string radius. Reciprocal torque, rotational spring storage and dashpot loss
therefore enter the same solve as bending, all other strings and the shared
body/cavity. Its force limit bounds `|torque|/r`. Body acceleration supplies the
existing pressure observer after bow release; twist receives no direct audio
gain. The new binding must precede playback and string motion. Original source,
body, gas, pressure and interface bindings survive multiple-bow scheduling.
G2 checks an independent two-oscillator angular-spring solution; G3 exercises
torque-only release, shear and body-density substitutions, original spatial
footprints and joint bowing; G0 covers rotational preload energy and refusals.
This does not establish nonlinear bridge contact, moving endpoints, local
plate-slope convergence, large-angle torsion, exterior radiation loading, or
measured violin/cello/bass fidelity.

### Material string bank independently directed translational bridges

`bowed_string::runtime::reciprocal::bank::polarization::StringBodyTranslationConnection`
and `ReciprocalStringBank::with_translational_body_connections` add independent
centreline translation laws to the existing material string bank. Each law
retains an original string index, interior station, signed unit
`TransverseDirection`, original-body `ModalAttachment`, nonnegative stiffness
and viscous resistance, and signed rest extension. The body map must measure
displacement along the SAME physical signed direction as the string map.
Body indices address the supplied bodies before strings are prefixed.

The relative coordinate is `direction dot y_string - x_body - rest`.
The original material source supplies the string projection with its actual
gas-loaded transverse mass. Its transpose supplies reciprocal force in the
same `CoupledModalSystem`; spring storage and dashpot loss enter the existing
work ledger. Original scalar bridges remain active, so distinct directions
and body maps can load both transverse planes without duplicating the
material string or prescribing modal detuning. Centreline translation has
zero direct torsional projection. Neither a complete rigid bridge nor its
contact law is inferred.

The retained source exposes `translational_body_connections()`. Cold gas,
cavity and rotational-bridge rebuilds retain those declarations and lower
their string maps against the current sourced realization. Original body
maps, material memory, gas state, exterior radiation, pressure observation
and dry bow interfaces stay with their existing owners. Binding accepts
one complete nonempty set before motion or playback and before nonlinear
storage; an empty set is a no-op. A one-plane source refuses any nonzero
second direction component. Scalar omitted-mode residual bridge flexibility
refuses this composition. Invalid stations, body indices, nonfinite or zero
body maps, nonpassive laws and owner budget violations refuse without
publishing a replacement system.

The `material_strings` frontend accepts an optional
`translational-bridges COUNT` block after thermal-support and rotational-bridge
records, before cylinder gas loading and `cavity`:

`translational-bridge STRING STATION SX SY SHELL_COMPONENT NODE GX GY GZ K C REST`

It requires an original material shell or assembled-shell body. The adapter
dots a global unit direction with the selected physical node's XYZ
translation in every retained `physical_modes` vector, preserving the original
mass normalization and eigenvector signs. Assembly components share ONE
assembled modal basis; their physical component index selects geometry,
while the modal attachment uses body zero. Recompiling after a body-material
substitution regenerates this projection from the new material solve.

String-local and body-global directions declare an orientation mapping.
The input has no string-to-world frame, so consistency between those
directions, including across connections, remains the geometry author's
responsibility. Unit-vector admission does not verify that mapping. The
frontend admits at most 16 additional rows and requires a retained second
plane whenever `SY` is nonzero.

Core regressions cover transverse-basis covariance, reciprocal second-plane
pressure, material substitution, zero direct torque, energy closure, cold
admission and gas/rotational attachment order. Frontend regressions project
the actual shell node, retain gas/cavity bindings, compare callback sizes
and change the body material. `directional-bridge.material-strings` is a
synthetic integration fixture, not measured instrument or bridge
qualification. Full three-dimensional unilateral bridge contact, moving
geometry and measured violin/cello/bass fidelity remain separate requirements.

### Omitted bending flexibility at fixed material bow contacts

Fixed bows may select `bow-residual-flexibility exact-static-green` after
the matching bridge `residual-flexibility` record. The material string bank
derives omitted bending compliance from each admitted string's tension,
bending stiffness, span and retained spectrum, reusing the bridge/point-force
Schur complement. Reciprocal string/body maps and all cross-bow terms enter
the existing simultaneous friction solve. Tail storage is `F^T G F / 2`;
`delta_F^T G delta_F / 2` is separately reported backward-Euler numerical loss,
not material or friction dissipation. Traction history publishes only with
the complete accepted step. This quasi-static approximation requires modal/time
refinement and does not restore omitted resonances or torsional modes. Fixed
point or uniform-ribbon bows on linear banks with one simple bridge per string
are supported. Source-defined shear films compose with this same residual,
including elastic/plastic/Maxwell history and series compliant drives. Each
film receives the actual interval slip `u dt - b^T delta_q - G delta_F`, further
reduced by any series-drive displacement. Omitted-mode storage and its
backward-Euler numerical loss remain separate from film storage and physical
heat. Finite-cell thermal trials retain the same accepted film and traction
histories, then publish the mechanical and thermal endpoints together. Either
attachment order reaches the same solver.

Fixed linear resolved fibers also compose with this string residual. Each
contact retains its original receiving fiber and adds the condensed bridge/body
correction to the same signed virtual-work column; the resulting column measures
slip and distributes reciprocal forces. Only omitted string bending modes enter
the residual compliance. Retained fiber modes, material losses, finite film heat
and the original body stay with their existing owners. String truncation, fiber
span/width discretization and time step require independent refinement. Tests
compare an independent five-oscillator solution, an explicit omitted-mode sine
series, receiving-fiber and body reactions, and fiber-modulus/film-heat-capacity
substitutions with thermal rollback.

Normal contacts, moving maps, nonlinear fiber strain, extra bridge maps,
nonlinear string storage, Mindlin and dry-interface source histories still
refuse. This composition adds no film opening/transport or measured rosin law.
Focused tests compare independent two-port elastic/plastic/Maxwell and compliant
drive solutions, retained pressure/history under reversal and rejected opening,
and heat-capacity substitutions without heating numerical tail loss. A native
material-string input reaches the same coupled thermal PCM path with identical
output across callback sizes. These are manufactured integration cases, not
measured instrument qualification. The optional Cannone experiment preserves
physical data and existing acoustic acceptance gates.

### `string_specimen`

`acoustic_realize::LinearMaterialStringRuntime` incrementally advances an explicitly
linearized pinned string from a `ResolvedStringSpecimen`, using the existing
exact-ZOH modes, loss admission and compact strip acceleration/jerk observer.
Material/operator state, modal displacement/velocity, observer history and sample
position publish at one epoch after all admission checks and a final `Cx` poll.
`step` returns a `MaterialStringFrame` pairing that epoch and the acoustic frame
with `StringDissipation`: material Kelvin-Voigt heat, oscillating-cylinder air
loss, and authored damping without an identified thermal destination. Coefficients
are retained when the modes are constructed and rebuilt with each rebind.
For constant additive viscous coefficients during a sample, each mechanism's
energy is its coefficient fraction times the mode-local viscous loss: all share
the same integral of squared velocity. This follows the viscous work law in
[MIT 8.01, section 23.5](https://ocw.mit.edu/courses/8-01sc-classical-mechanics-fall-2016/mit8_01scs22_chapter23.pdf).
The existing Rayleigh option still replaces the entire damping model; its loss
is authored, with no silently added air or material contribution. Tiny negative
mode-local differences and undamped drift remain in a signed roundoff residual,
not in heat. That residual must fit the combined modal/aggregate/reduction
roundoff allowance before the runtime commits. Diagnostic microphone pressure
adds no second radiation-energy debit. The returned partition does not itself
debit reservoirs or establish a closed thermal feedback/entropy model.
Refused steps and rebindings preserve the complete accepted runtime. Cloning an
accepted runtime retains the acceleration history needed for identical continuation.

`project_point_forces` maps prescribed transverse loads in newtons at material
stations to the current mass-normalized sine basis. `point_motion` uses the same
participation to return accepted displacement [m] and velocity [m/s]. The transpose
pair preserves power and held-load work; endpoint participation is exactly zero
for pinned supports. The resulting forces feed the existing mechanical and thermal
steps, including after refinement and material updates. Projections must be used
before the next epoch change. This does not supply contact laws, opposite-body
reactions, moving-load work or compliant-support dynamics. G1 checks an independent
physical oscillator and work; G3 checks multipoint power/work through thermal
updates; G4 checks endpoints, invalid inputs, overflow, cancellation and retry.

`initialize_point_load_equilibrium` connects the same physical force projection
to the modal owner's static solve `q_k = f_k / omega_k²`. It initializes only an
unexcited runtime before its first sample, establishes zero prior acceleration,
and advances the epoch without changing material state or the sample clock.
Holding those forces retains equilibrium; releasing them excites the same live
mechanics and pressure observer. Stored initial energy assumes a settled loading
history before the simulation window; it is not a simulated actuator-work or
contact receipt. Existing motion cannot be overwritten. G1 checks static compliance
and release response; G4 checks repeat, invalid-force and cancellation refusals.

`rebind_ambient` resolves a new prescribed gas through the existing gas owner
and rebuilds air drag and compact observer weights atomically. Material receipts,
modal coordinates, previous acceleration and accepted sample count are retained;
the epoch advances with zero parameter work because this model has no gas-added
mass or gas-dependent elastic storage. The thermal string exposes the same update
with expected-epoch checking and unchanged material enthalpy. Boundary heat must
still enter through thermal transport. This is a piecewise prescribed ambient
change, not an atmosphere transient, fluid-loading work or propagation-delay solve.
G1/G3 checks use cylinder drag and compact pressure formulas including retained
jerk history; G4 checks invalid gas, stale epochs and cancellation preserve retry.

Rebinding requires the exact same immutable card, length, linear density and
sine-mode count. This makes the mass-normalized basis identical, with no sorted-mode
matching, displacement reset or sound-level rescaling. Tension, EI, damping and
observer weights rebuild together; the source state remains inspectable by other
consumers. Signed parameter work is the change in modal elastic energy at fixed
displacement/momentum, checked against an explicit caller budget. It is returned
for the thermal/external parameter owner's opposite transfer, not booked as damping.
The axial prestress reservoir, heat transport and closed thermomechanical balance
are not supplied here. Coefficients are externally scheduled and piecewise constant.
The acoustic observer retains the preceding accepted acceleration across updates;
abrupt coefficient changes can produce a transient, without a crossfade.

`substitute_material` additionally accepts a different source card and line mass
at exactly the same length, circular radius and pinned sine-mode count. It uses
the same coefficient/loss/observer rebuild, with an explicit canonical transfer:
`q_new = sqrt(mu_new / mu_old) q_old` and
`p_new = sqrt(mu_old / mu_new) p_old`. Physical displacement and transverse
momentum remain continuous; velocity generally changes. The reported external
parameter work includes both kinetic and elastic energy changes and must fit the
caller's budget. Previous acceleration is rescaled into the new mass normalization
so the pressure observer retains its physical history. Physical forces must be
projected again after the epoch changes. The compact observer includes endpoint
acceleration changes but does not reconstruct distributional radiation from the
instantaneous velocity jump; that needs a time-resolved transition.
G1 checks an independent forced oscillator,
canonical momentum and pressure including the even-mode jerk term; G3 checks a
substitution round trip and identical-source continuation; G4 checks complete
rollback after geometry/basis, work-budget and cancellation refusals. This is a
prescribed material intervention, not chemical conversion, melting, mass transport
or heat deposited in either material. The thermal wrapper retains its stricter
same-card binding and does not expose this substitution.

Explicit `refine_modes` increases the retained sine basis at fixed geometry without
advancing time. Existing coordinates and observer acceleration history are retained;
added modes start at rest with zero previous acceleration. Shared admission enforces
mode, Nyquist and material-frequency limits before publication. The thermal wrapper
checks the expected epoch and preserves enthalpy. G1/G3 checks continuation energy,
pressure and independent forced response; G4 checks atomic refusals and cancellation.
This does not reconstruct omitted history or estimate truncation error, and does not
support coarsening or transfer between different geometries or eigenbases.

The rung omits transverse-motion-induced axial stretching, moving supports, extra
polarizations, Prony memory, basis/geometry transfer, phase change and radiation
backreaction. It uses prescribed piecewise ambient gas and fixed listener distance, with no atmospheric
path filter or common propagation delay. G1 uses a synthetic thermal-expansion card
and an independent damped-oscillator/pressure calculation. G3 cycles a same-card
density/temperature update at prescribed mass through changed cross-section, bending
stiffness and pressure without changing the retained sine basis. G4 exercises material-work
and late observer refusals, cancellation and exact suffix replay. This is focused
implementation evidence, not measured string fidelity or full MR05/MR34 completion.
The loss-partition G1 test independently integrates squared analytic modal
velocity and feeds only material heat into the existing equilibrium enthalpy
owner. Its synthetic constant-coefficient solid absorber checks the combined
mechanical/thermal/air balance and rejects duplicate heat, air credited to the
solid, and a missing mechanical term. It stays wholly on the solid branch;
thermoelastic feedback and phase dynamics are not tested. G3 changes loss laws
without resetting vibration and checks that the returned destinations change.
Existing G4 refusal/retry tests compare the complete frame, including loss fields.

`acoustic_realize::ThermalMaterialStringRuntime` owns that incremental string,
its immutable material card and a uniform equilibrium enthalpy state. It closes
the dissipative temperature/viscosity feedback. Its ordinary constructor admits
the prescribed-radius or prescribed-mass, fixed-tension, constant-density,
constant-Young's-modulus solid rung. Those
restrictions keep elastic and axial prestress storage independent of temperature;
there is no omitted stiffness-change or expansion work. Density must agree with
the solid branch of the thermal chart. The initial absolute temperature and card
must agree across both owners. Original query-axis schemas and selected claims
are retained; subsequent queries pin those claims, with no automatic substitution.
Constant properties and linearly interpolated temperature-only laws are admitted;
exact-only temperature tables do not authorize an evolving path between samples.
Viscosity must be nonnegative throughout its supplied source curve, including
interior knots whose values cannot be inferred from positive endpoints.
For `FixedMass`, the prescribed mass and `m/L` remain exact through thermal
updates; constant density and loaded length keep the resolved radius and basis
unchanged. Thermal area and volume use `ResolvedStringSpecimen::radius_m()` and
its cross-section, rather than reconstructing mass through a rounded radius.
This ordinary constructor does not admit changing density, expansion or geometric
state transfer. The source-bound thermoelastic constructor below also keeps mass,
section and loaded length fixed while evolving its thermal eigenstrain.

`try_new_thermoelastic` admits the original sourced `FixedThermalExtension`, a
scalar expansion coefficient and scalar or globally affine temperature-dependent
Young's modulus. Every original modulus knot must belong to that same affine law;
derivative jumps require a different caloric model. The original material and
expansion owners resolve the complete temperature path, retaining source claims,
reference temperature and reference stress-free length. At each accepted sample,
thermal expansion changes the actual tensile force and modal operator without
resetting displacement, velocity, observer history or the accepted sample count.
Slack, changing density/geometry, melting and unsupported source paths refuse.

The axial Helmholtz storage is `Psi_a = A L_ref E(T) e(T)^2 / 2`, with
`e = (L-L_ref)/L_ref - integral(alpha dT)`. The retained native transverse
storage is `Psi_v = sum(omega_n(T)^2 q_n^2 / 2)`. Mechanical entropy is
`s = -partial_T(Psi_a + Psi_v)` at fixed supports and modal coordinates, and
mechanical internal energy is `K + Psi_a + Psi_v + T*s`. The public
`mechanical_internal_energy_j` reports this same storage for coupled consumers.
Positive total fixed-coordinate heat capacity is checked on every original
caloric segment, including the amplitude-dependent modal term; the exact
quadratic extrema of `T*Psi_TT` are checked for the admitted affine law.

Mechanics advances once at the accepted temperature. At its resulting modal
coordinates the existing bounded scalar feedback solver closes reversible heat
`Q_rev = -delta(Psi_a) - modal_parameter_work - delta(T*s)` through the original
caloric/transport owner. The heat includes entropy change due to modal motion,
not just the endpoint temperature change. Thermal callbacks are repeatable
`FnMut` proposals from the same accepted state; only the final report/reservoir
candidate survives. The ordinary constant-stiffness path still calls transport
once. The final frame exposes axial, modal-parameter and entropy-energy changes
separately, updates the actual modal energies and pressure, and checks the full
first law. Material dissipation remains distinct from signed reversible heat;
neither numerical loss nor diagnostic acoustic pressure is deposited as heat.
The electrical and finite-gas paths share this owner, including reciprocal
conduction at the actual solved solid and gas endpoint temperatures.

This remains a first-order partitioned small-strain model with a fixed section,
uniform temperature and diagnostic compact radiation. Energy closure does not
prove exact entropy conservation, timestep convergence or instrument fidelity.
`tests/string_specimen/thermoelastic.rs` defines synthetic live-pitch, independent
internal-energy, electrical-source and refusal/retry checks. Native compilation
and execution were unavailable in the implementation environment; syntax and
energy-algebra checks are not reported as native test passes.

Each sample advances mechanics at its initial viscosity, credits material loss
and explicit external heat to total enthalpy, then resolves the new material
temperature and rebuilds damping. Final pressure and retained observer history
use the updated acceleration. This is a first-order partitioned coupling; exact
mechanical substeps do not make the coupled solution exact. A positive caller
temperature-increment budget constrains every accepted step, and temporal
refinement must assess coupling error. Air loss and mechanical roundoff remain
outside material heat. The returned first-law residual includes work, both storage
changes, outgoing losses and numerical loss, with a scale-aware arithmetic
allowance including the enthalpy reference. It is not a physical-error certificate.
All candidates publish at one epoch after the final cancellation poll. Reusing
an accepted input epoch refuses before consuming another heat input or sample;
thermal/material refusal after the mechanical substep leaves both owners intact.
Cloning preserves exact continuation, including thermal state and source pins.
The direct `step` takes external heat in joules. `step_with_current` instead
resolves the existing `fs-material::conductor` from the specimen's pinned
electrical-resistivity claim, current temperature and actual length/area. An
ideal source deposits `I² R(T_initial) dt` once in the external heat account;
material vibration loss remains separate. The returned `ElectrothermalStringFrame`
retains initial/final sourced resistance and the same accepted thermal/mechanical
epoch. Final electrical resolution precedes publication; any refusal preserves
the accepted state and retry position. This is first-order partitioned uniform
Ohmic heating with insulated boundaries, not circuit dynamics, skin/contact
resistance, a charge-field solve, Lorentz forces or motional emf. Synthetic
G1/G3/G4 tests check analytical temperature refinement, current-sign symmetry,
resistance/damping/pressure response and atomic retry, not experimental fidelity.
`step_with_current_and_thermal_transport` composes the same current source
with the existing transport callback. The thermal owner receives vibration
plus Joule heat as internal deposition and returns actual boundary heat.
`ElectrothermalStringFrame<R>` retains Joule heat, boundary heat and the
concrete transport report separately; its coupled external-energy account is
their sum, so electrical work is neither lost nor counted twice. Both paths
share current admission, and final electrical resolution still precedes
publication. Zero current reduces to the ordinary thermal transport trajectory.
The independently prescribed fluid/enclosure reservoirs retain the existing
uniform-body, Biot and infinite-reservoir restrictions.
`step_with_voltage` and `step_with_voltage_and_thermal_transport` share these
same transactions under `OhmicDrive::Voltage`. Each sample freezes the current
`V/R(T_initial)` and deposits its Joule energy; the next sample uses the updated
resistance. The frame retains the authored drive, initial sample current and
current demanded at the final accepted resistance. Final electrical evaluation
must succeed before publication. Positive temperature coefficients therefore
reduce fixed-voltage heating while increasing fixed-current heating. The
partition remains first order; no within-sample circuit solve is claimed.
`step_with_thermal_transport`
accepts a thermal proposal through the shared `fs-material::phase` input/result
carriers. Passing `|input| environment.advance(cx, input)` with
`fs-conduction::lumped::LumpedThermalEnvironment` uses the existing
`solve_lumped_enthalpy` for one
backward-Euler step over the same audio sample. It derives the exposed cylinder
area `2 pi r (L+r)`, volume and Biot length `V/A` from the resolved specimen,
with no independently authored thermal mass or geometry. Convection and surface
radiation act on this same complete boundary; supports are thermally insulating.
The environment is a prescribed infinite, uniform reservoir and enclosing
radiation field with independent fluid and enclosure temperatures. Opposing
convection and radiation fluxes contribute their signed heat once. It neither changes nor debits the separately
declared acoustic gas. Transport may be declared or card-backed; the conduction
owner enforces its existing card identity and temperature-coverage rules.
The callback must not mutate or debit an external source before acceptance.
`step_with_thermal_state` also accepts caller-owned state by exclusive borrow.
Its callback reads that accepted state and returns an owned candidate alongside
the normal thermal proposal. Only after all coupled checks succeed does the
method replace the caller's state. A finite support therefore advances with
the string without a separate caller commit. G3/G4 finite-support integration
checks heat debit, temperature/damping/pressure feedback, late foreign-chart
refusal and exact retry of both owners. The transport owner validates its
opaque state; callbacks must not mutate shared interior state or external
resources. This is synchronous publication through exclusive borrows, not a
cross-thread or durable transaction.
`step_with_drive_and_thermal_state` provides the same publication boundary
for either `OhmicDrive::Current` or `OhmicDrive::Voltage`, reusing the existing
electrothermal path. External thermal state remains a candidate until final
resistance/current resolution also succeeds. The callback receives Joule plus
material heat; its boundary transfer excludes both. Finite-support tests check
this separation, current/voltage callback parity, cancellation and retry.
The coupling owner verifies that the proposal uses its original phase chart
and finite heat/residual budgets, then checks the coupled energy balance. This
keeps fs-conduction out of fs-couple's production dependencies: FEEC's optional
topology feature already depends on fs-couple. The cross-crate test has a
development-only conduction dependency.

Material loss enters once as mean internal power for that sample. Boundary heat
uses the actual endpoint convection/radiation powers, not a subtraction from
the computed enthalpy change that could hide thermal solve error. The returned
transport trajectory retains those powers, its Biot value and energy residual.
The caller's thermal energy-residual budget is checked separately and appears
as `thermal_solve_tolerance_j`, zero on the direct-heat path. Coupled first-law
acceptance uses that budget plus the reported arithmetic allowance; neither
certifies time-discretization error. A late thermal, phase or material refusal
preserves mechanics, temperature, source receipts, epoch and observer history.

The fixed corpus Biot ceiling is applied, including the conservative radiative
coefficient over the chart. This is a necessary uniform-body check, not proof
of initially uniform temperature or homogenized internal heating. Rapid local
sources and nonuniform support transfer require spatial transport regardless
of step size. The analytical convection reference and lumped assumptions are
described in [MIT's transient heat-transfer notes](https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node129.html).
The thermal chart retains its existing source authority. Finite-support contact
uses the conduction adapter's fixed resistance and uniform-body limits; spatial
conduction, entropy production certification, thermal
expansion, changing stiffness/density, and solid-liquid transfer are not claimed.
G1 compares the coupled trajectory with an independently integrated continuous
ODE and checks temporal convergence. Synthetic heat capacity and viscosity make
feedback observable in a bounded test; they are not measured material data.
G1 additionally checks analytical heating/cooling and first-order refinement
for two radii, plus separate and combined convection/radiation inputs changing
the material damping and diagnostic pressure. G4 tests Biot, thermal residual,
source-domain, phase and cancellation refusals, duplicate epochs, retry and
cloned continuation. These exercise the existing solver directly, not a mock.

`with_uniform_circular_material_state` lowers one homogeneous isotropic circular
specimen to the existing `PrestressedString` realizer. It consumes the shared
`ResolvedMaterialStatePoint` API, requiring only positive SI density and Young's
modulus. Area, area second moment, linear density, mass, EA, EI, and diameter
follow from the same radius and length. All source receipts and the complete
material state remain available in `ResolvedStringSpecimen`; its specimen hash
binds geometry and material resolution, excluding mechanical loading and solver
options. It is not the complete scenario identity.

The original binding holds the template's geometry and tension fixed.
`with_uniform_circular_material_and_prestress` additionally takes an explicit
`StringPrestress`: fixed tension, fixed extension from a stress-free length, or
target primary undamped fundamental frequency. Fixed extension computes `T=EA*strain`
and enforces the caller's declared linear strain limit; that assumption is not a
measured yield boundary. Target frequency inverts the existing pinned-pinned
Euler–Bernoulli dispersion relation, including EI. Slack/compressive solutions and
extension/target prescriptions on the moving-end basis refuse. The prescription
remains separate from the material receipts, available through `prestress()`.

`with_uniform_circular_material_and_constraints` also accepts a
`StringGeometryConstraint`: fixed radius or fixed total mass at the template's
current loaded length. Fixed mass derives `r=sqrt(m/(rho L pi))`, retains the
prescribed mass, and sets `mu=m/L` directly. It never reconstructs either mass
quantity through the rounded radius. Thus changing density at fixed mass/length
keeps exactly the same numerical modal basis while changing EA, EI and observer
width. The radius-derived area and `rho A` agree with these mass constraints to
floating-point roundoff. Geometry and prestress constraints are independent and both remain
available on the result, separate from material authority. The specimen hash
binds resolved geometry, mass and linear density, not the authored comparison
policy. The v2 identity distinguishes prescribed masses even when rounding maps
them to the same radius.

`with_uniform_circular_thermal_extension` consumes the existing
`IntegratedIsotropicThermalExpansion` and current elastic resolution from the
same exact card, query coordinates and axis conventions. It computes the
small-strain eigenstrain relation `T = EA * ((L-L_ref)/L_ref - integral(alpha dT))`,
using one reference configuration, as in
[MIT 16.20, section 3.7](https://web.mit.edu/16.20/homepage/3_Constitutive/Constitutive_files/module_3_no_solutions.pdf).
The current geometry/mass constraint remains explicit; both elastic and expansion
receipts survive in the result. `FixedThermalExtension` records the mechanical
prescription, while `thermal_expansion()` distinguishes the sourced path from
raw authored eigenstrain. A subsequent ordinary binding requires the thermal
path to be supplied again; it cannot silently retain old thermal evidence.
Total, thermal and positive elastic strains must stay within the authored limit.
Slack, compression and domain exit refuse; no buckling solver is implied.
Fixed tension remains a separate force-controlled prescription. G1/G3 check a
piecewise-linear expansion curve across temperature knots, heating/cooling force
changes, mass and receipt retention, and pitch from the real pressure solver.
G0 checks current-state/source mismatch, invalid reference length, unsupported
moving supports and strain-domain refusals. This initializes independent uniform
specimens; it does not evolve temperature, exchange heat, transport strain or
solve lateral contraction during playback.

All policies preserve the selected loss law, supports, and modal budget. These compare
specimens at resolved states; thermal evolution, transverse contraction,
yielding, anisotropy, or melting are not solved. Target pitch describes an undamped
linear eigenfrequency, not a damped or finite-amplitude pressure frequency. Diameter
supplies the existing compact strip-radiator width, without a new claim of accurate
circular-wire radiation. The density and modulus requirements use exact
dimension-only quantity schemas; equal-dimension semantic kinds or waveform
forms cannot be erased when lowering to the solver. Admission refuses missing
or quantity-schema-mismatched properties,
invalid mechanical inputs, and nonfinite/zero derived products before publication.
The finite scalar calculation is deterministic; no inner-loop database access is
added. `tests/string_specimen.rs` checks G0 admission, G1 mass,
EA/EI and beam modes, and G3 geometric/material scaling plus pressure-measured pitch
through the actual acoustic realizer. G1/G3 cases also cover extension-derived
tension and material changes under explicit extension/target-frequency constraints,
plus fixed-mass diameter/stiffness changes and their bending-sensitive pressure pitch.
Synthetic data check numerical behavior, not physical validation of a measured material.

`ResolvedStringSpecimen::with_kelvin_voigt_bending_loss` selects a constant
Kelvin–Voigt bending law using the retained `kelvin_voigt_bending_viscosity`
property [Pa s]. Density, Young's modulus and viscosity must be dimension-only
validity-wide scalar constants or temperature-only curves at fixed positive
absolute T. Frequency samples are not promoted to a broadband law. The
viscosity must declare a frequency domain, intersected with any
frequency restrictions on density and modulus. Other state coordinates stay at
their admitted query point. The descriptor carries `eta I`, the intersected band,
and the material bundle identity. Rebinding an already viscous specimen recomputes
all three from the new geometry and material. No source uncertainty is improved.

All three material loss binders (Kelvin–Voigt string/plate and Prony string)
use `ResolvedScalarProperty::angular_frequency_band_rad_s`. One static typed
frequency coordinate can have any name and be expressed in Hz or rad/s;
legacy untyped `omega` remains rad/s. The normalized band intersects every
participating coefficient's applicability before the existing modal checks.
Multiple frequency coordinates refuse because they may be independent physical
variables. No source receipts or interpolation decisions are rewritten. G3
tests compare actual pressure bit-for-bit across the two unit conventions and
legacy inputs; narrowing only density support rejects an out-of-band string.

The shared linear/Kirchhoff–Carrier/moving-end modal path projects
`mu y_tt - T y_xx + EI y_xxxx + eta I y_xxxxt = 0`, giving
`zeta_k = eta I k^4/(2 mu omega_k)`. Thus tension contributes stiffness but no
viscous loss. This is the constant-coefficient Kelvin–Voigt term in
[Sakthivel et al. (2023)](https://arxiv.org/abs/2301.07931), with tension added
and projected onto the existing sine/cosine bases. The explicit law replaces
the authored constant modal ratio; simultaneous Rayleigh
or nonzero authored damping refuses. Air resistance is still added separately.
Every retained reference frequency, including the secondary polarization, must
lie in the descriptor band. G1 checks `eta I` and rebinding; G3 checks actual
pressure decay against an independent diameter-form expression in all three
solver paths, including zero-viscosity air-only decay. G0 covers unavailable
data, wrong quantities, sampled frequency curves, duplicate losses and out-of-band modes.

This Kelvin–Voigt choice models linear bending strain-rate loss only. It does not supply axial
viscosity, material memory, amplitude-dependent loss, time-evolving material
state, thermoelastic/joint/radiation dissipation, or measured material validation.
The modal applicability check concerns reference frequencies; nonlinear
frequency shifts and generated harmonics are not an in-band certificate. The
existing moving-end cosine approximation and compact observer remain unchanged.

`with_prony_bending_loss` binds a causal generalized Maxwell bending spectrum
from explicit modulus [Pa] / relaxation-time [s] property pairs and
`equilibrium_young_modulus` [Pa]. `with_standard_linear_solid_bending_loss`
selects the single `relaxing_bending_modulus` / `bending_relaxation_time` pair
through that same binding. The explicit equilibrium modulus must equal the
Young's modulus used to derive EA/EI; a sampled dynamic modulus cannot silently
be treated as an equilibrium value. All selected coefficients, including density,
must be dimension-only validity-wide constants or temperature-only curves at
fixed positive absolute T. Every relaxation-time claim must
declare a frequency band, intersected with all other coefficient bands. An empty
spectrum is explicitly elastic and requires an equilibrium-modulus band. Geometry
or material rebinding resolves every retained source pair again, recomputing
each delta EI and replacing the source identity and time constants. Missing or
mismatched selectors refuse rebinding; duplicate source pairs refuse binding.
These are resolved inputs, not measured-material validation.

Each mass-normalized modal coordinate receives generic `fs-phs` relaxation arms
with `a_j = delta EI_j k^4/mu`: `q'' + c q' + omega0^2 q + sum_j a_j(q-v_j)=force`,
`v_j'=(q-v_j)/tau_j`. This is the internal-strain generalized Maxwell formulation
described in [Bleyer's viscoelasticity tutorial](https://bleyerj.github.io/comet-fenicsx/tours/nonlinear_problems/linear_viscoelasticity_jax/linear_viscoelasticity_jax.html),
projected onto the existing beam modes. Nonlinear axial storage, contact ports,
body joins and pressure observation use the same augmented pHS state. A held
pluck initializes every `v_j=q`, representing full relaxation before release.
Positive-strength branches add memory states in branch-major, mode-major order;
zero-strength branches retain source selection but add no state or active pole.
The owned spectrum is borrowed throughout the audio loop. These
material arms produce storage dispersion and dissipation together; no PCM decay
envelope or substitute modal damping ratio implements it.

Admission checks both relaxed and instantaneous reference frequencies against
the material band and Nyquist guard, using the sum of all branch stiffnesses
and including fixed-end detuned polarizations.
Moving-end secondary polarization refuses on this path. Kelvin–Voigt, Rayleigh
and nonzero authored internal loss cannot be combined with the Prony descriptor;
air drag remains separate at the relaxed reference frequency. G0/G1 tests cover
sample-clock admission (`dt/tau_j <= 0.1` for every nonzero arm and
`dt*omega_instant <= 0.1`, otherwise increase sample rate),
source/quantity/band refusals, geometry rebinding and actual pressure pitch/decay
against independent characteristic-polynomial roots in fixed linear, nonlinear
and moving-end realizations. Multi-arm G1 pressure is compared with independently
integrated hereditary-strain dynamics and checked for second-order refinement;
G3 checks equivalent branch splitting/permutation, zero-arm omission, detuned
polarizations, and body/contact composition in fixed and moving-end paths.
This remains linear bending memory at a fixed state. No measured spectrum fit,
axial creep, evolving temperature, phase change, arbitrary pre-release loading,
or in-band certificate for nonlinear harmonics is supplied.

### `driving_point::regional::field`

`RegionalDuct::with_exterior_gas` now binds side-hole radiation as well as the
terminal load and observer. Every regional run retains its own chimney gas
and receives the same explicit exterior through `fs_duct::exterior`; pressure
and volume flow join the two media without a second end correction. Open and
vented holes may therefore contain a different mixture from the atmosphere.
This also reaches the continuous-field response/reflection, baked terminal
loads and scheduled reed renderer. A closed main outlet does not disable
radiation through open side holes. An omitted exterior retains the original
local-gas convention and identical media retain their original bits. Static
pressure mismatch still refuses, and the segment-gas and midpoint bindings keep
their varying-medium mutual-hole refusals; the continuous-field binding applies
the stated geometric law below. Focused G1/G3 tests cover physical side-hole transfers, upstream gas
regions, continuous fields and changes to both reed motion and microphone
pressure. This is frozen pressure-matched acoustics, not a mixing/mean-flow
model or experimental validation of a complete instrument.

`RegionalDuct::from_gas_field` discretizes the shared `fs_material::gas_field::AxialGasField`
at explicitly refinable midpoint cells (at most 64 total), retaining the
existing cylinder/cone TMM and pressure/volume-flow interconnection. It adds no
propagation or loss law. The exact x=0 state feeds source density and the inlet
characteristic impedance, and x=L feeds the existing terminal radiation load;
an independent exterior can still be supplied explicitly. Uniform fields keep
the original complete geometry and arithmetic. Nonuniform tone holes and smooth
radius profiles refuse pending the corresponding continuous-field duct support.
`ReedBoreVoice::new_with_gas_field` feeds this same regional line into the existing
reed runtime and compact-jet observer. No instrument constants, tuning law,
evolving breath or calibrated exterior observation are added.

The focused field tests compare uniform responses bitwise, exercise cone
radiation and endpoint states, and refine a supplied 1 K axial temperature
gradient over 4/8/16/32/64 cells. Fundamental resonance convergence is second
order; the third-to-first register interval is checked against independent
energy-weighted perturbation. A rendered-pressure test checks uniform legacy
identity, exact source density and the gradient's nonzero audio effect. These
are numerical tests of prescribed fields, not measured instrument qualification;
native results are recorded separately.

`RegionalDuct::from_continuous_gas_field` uses the existing fourth-order
Magnus operator for simultaneous smooth geometry and primitive gas fields.
Its optional `WallPin` now reaches the same local wall-admittance owner as
uniform cylinder/cone ducts at every gas/geometry quadrature point. Responses,
characteristic-line feedback and baffled outlet pressure retain that common
wall binding. The wall's explicitly supplied per-area mass/stiffness/resistance
remain constant; actual radius and slope set the surface coupling. This is a
locally reacting wall, not an axial shell, evolving wall temperature or a
nonlocal wall-radiation model. The continuous binding now admits tone holes
at their exact local gas state through the original homogeneous-chimney and
junction laws. Closed holes and any number of open/vented holes may lie in a
gradient. Every open hole is built in its own local gas and retains every
mutual correction, Nederveen's geometric near-field series length, which does
not depend on the gas: no averaged medium is substituted and no term is omitted
(the term is at most about 0.06 mm, about 0.3 cent, on the clarinet family and
its gradient-dependent part about 0.002 cent; the law is stated in the
`GasFieldDuct::new` doc). The segment-gas and midpoint bindings keep their own
refusals. This also reaches the
existing scheduled reed renderer without a separate propagation path.
Open/vented chimney flows now cross this continuous-field binding into the
same bell-and-hole microphone realization. Axial delays use the original
gas-field quadrature and each chimney's local sound speed. Feedback, chimney
mouth loads and exterior microphones use one homogeneous exterior: the explicit
binding, or the actual outlet gas when omitted. The same choice applies to
regional interfaces, arbitrary terminal loads and characteristic-line shortcuts. Selecting a continuous field cannot
silently turn an open hole into a non-radiating source. Closed holes remain
loads only. This retains the existing sampled-transfer and infinite-baffle
approximations, rather than adding a new radiation or transport law.

`RegionalDuct::with_open_hole_interaction` and the optional reed-card record
`open_hole_interaction tmmi EPSILON` (immediately before `compile_limits`, on a
`gas_field continuous` binding with no molecular bore; EPSILON 2 as the paper
prints it) opt in to the Lefebvre, Scavone and Kergomard (2013) coupling of the
open holes' and the bell's radiation that `fs-duct` documents for
`GasFieldDuct::with_open_hole_interaction`: feedback, reflection and the
observed hole flows come from one coupled solve. It is off unless declared and
then every code path is the original. The distances are axial (an upper bound),
and a segment-gas or midpoint binding, a uniform gas field (which gives the
wrapper no continuous binding), a vent fraction, an ideal-open end or a baked
mouth load each refuse by name rather than being ignored. The receiver's
pre-arrival travel times and each hole's nonlinear orifice port keep the
independent chimney.

`RegionalDuct::with_printed_hole_lengths` and the optional reed-card record
`tonehole_lengths lefebvre-scavone-2012` (after `open_hole_interaction` if both
are given, before `compile_limits`; a `gas_field continuous` binding with no
molecular bore) opt in to the thick-wall hole lengths of Lefebvre and Scavone
(2012) that `fs-duct` documents for `GasFieldDuct::with_printed_hole_lengths`.
It is off unless declared, composes with the interaction, and refuses by name on
any binding without a continuous operator.

The optional reed-card record `constriction_noise SEED DUCT_RADIUS_M` (after the hole rows, before
`compile_limits`; both fields required, no default, once only) gives the reed channel the brass lip row's
duct-constriction noise (`turbulence_noise::ConstrictionNoise`; Oldham and Ukpoho 1990 as printed in Guedel 2005,
the trend line read from the figure): a seeded third-octave series source on the mouth side, of rms K(St) times the
held blowing pressure, with `St = f d / U_c`, `U_c = sqrt(2 P / rho)` (the Bernoulli speed of the total drop) and
`d = pi r / 2` of the card's duct radius. It is silent while the previous sample's mouthpiece flow was not forward,
and off unless declared: without it every sample is the original. It is an extrapolation in size, speed and shape,
from dampers and orifice plates in a 300 mm duct at 10 to 25 m/s to a reed slit; the printed law fixes neither the
jet speed nor the pressure drop for a viscous slit, and the readings available on the engine's slit span 40 dB.
A focused impulse regression covers the combined wall/field changes to
both incoming feedback and causal external microphone pressure.

### `acoustic_realize` / `pcm_wav`

`realize_assembly_with_gas_and_absorption` accepts a caller-resolved gas and
an explicit `GasAbsorptionModel`. Its observer always uses the selected law;
invalid ISO-air conditions refuse without a classical fallback. `RealizeOptions`
carries the same resolved gas and observer law alongside the bore discretization.
Ordinary assemblies, regional plate charts, `CompiledMaterialAssembly` and
`CompiledValveMaterialAssembly` all forward those options to one dispatcher.
Material compilation therefore preserves the caller's explicit absorption and
reed-bore resolution, including their admission failures. The gas is validated
once before stepping and supplies excitation, duct coefficients, body/cavity
loading and observed pressure; solid card temperature and prestress stay
independent. Missing absorption selection on a supplied gas chooses classical
Stokes-Kirchhoff by input semantics, never by coefficient equality with a preset
or an invalid humidity sentinel. The scenario-air entry point preserves its
legacy absorption arithmetic. `RayleighMedium::from`
copies density/sound-speed bits without reevaluating an atmosphere. G3 tests
cover explicit ISO/classical selection through source-bound uniform/regional
solid bodies and material valves, pressure/composition changes reaching actual
pressure histories, unchanged solid identities, retained bore refusals, invalid
ISO humidity, exact medium views and legacy history equality. Ideal-gas pressure
changes density/loading without an invented sound-speed or pitch retuning.
GasState itself carries no source receipts; arbitrary-mixture molecular
relaxation, changing composition during a render and experimental validation
remain outside these frozen-state paths.

File-driven bowed, reed and plate-valve performances expose
`from_bytes_with_gas_sources` and retain their authored `GasRecord`s. Ambient,
regional and exterior declarations use the same L3 reader. The CLI loads
bounded source artifacts itself: repeated `--gas-component SPECIES NASA9
TRANSPORT-CLAIMS`, explicit `--gas-water none|ID` and `--gas-phase
frozen|liquid-water`. A primitive reed render additionally supplies
`--gas-record FILE`; wind/ensemble gas records live in their performance files.
No thermochemical route falls back to a preset. Source-aware sidecars retain
composition, temperature, pressure and all three pack identities per species.
Omitting the new options preserves legacy numerical and sidecar paths.

The native `percussion` example uses this same bounded source loader for
drum, cymbal-profile, explicit shell-mesh and paired hi-hat `ambient` rows,
including both shell files referenced by a pair and `--gas-record` overrides.
The selected `Ambient` owns the original `GasRecord`; its exact admitted
density and sound speed supply exterior radiation and distributed cavities,
and `rho*c^2` supplies the existing drum volume spring. Composition, phase
admission and species/NASA-9/transport identities remain with the run and its
ambient receipt. The same registry resolves a hi-hat squeeze film's separate
`gas` row; source viscosity and p/T/R enter the original compressible film
owner, which retains the complete record and its existing isothermal storage.
No new pressure force, damping term, time loop or solid retuning is added.
Missing, ambiguous and out-of-domain sources refuse, while legacy declarations
preserve their numerical state. Focused native tests cover actual source bytes
and nested file routes, NASA-9-dependent drum pressure feedback with exact
replay, and source species mass in compressed real-skin gas pockets. This
extends M0.1-W1 input-to-consumer coverage; it does not qualify a full-band
percussion instrument or add evolving gas temperature/composition.

The same native percussion route accepts `head_material` in a drum specification:
side, thickness, installed tension, explicit loss token and the shared
`MaterialRecord` (`PACK STATE [AXIS=VALUE ...]`). The existing native 8 MiB card
loader is shared with plate inputs; each named pack is loaded once per drum.
`material_source::plate_section` resolves isotropic E, nu and density at the
same physical query before the original `TensionedDisk` assembly. The full pack,
source record and resolved state remain with the executing drum. Numerical
coefficient overwrites and missing, ambiguous, foreign-state or out-of-domain
properties refuse. Only static properties at fixed external conditions are
admitted; a frequency/strain/rate sample cannot masquerade as an elastic law.
Geometry and installed tension remain explicit, with no inferred thermal
expansion or retuning. Existing head relaxation may compose only through its
unchanged equilibrium-modulus and loss-ownership checks. Focused regressions
cover exact numeric-source parity, independent mass/bending scaling, actual
contact/nonlinear-head/cavity response and measured-pack file loading. This
extends M0.5-W1 without claiming measured drum-film damping or instrument
qualification.

`realize_assembly` turns an `fs_scenario::AcousticAssembly` into observer
pressure by composing nameless primitives: a prestressed Euler–Bernoulli
waveguide (`fs-nlmodal::prestressed_beam_omega` / Kirchhoff–Carrier), a
`bernoulli_aperture`, a `stribeck_friction` contact, a
`traveling_wave_line`, `ModalAcousticTime`, and compact modal radiators.
There is no instrument crate and no instrument algebra. A guitar or
clarinet is one filling of those objects.

- Ambient state comes from `GasState::try_new_moist_air(T, p, RH)` once,
  before realization. Excitation, duct impedance/resonance, radiation loads,
  and observer absorption share that state. RH=0 uses the original dry-air
  constructor bit for bit; RH>0 requires 253.15–323.15 K and water mole
  fraction at most 0.15, with model-domain failures returned as `Ambient`
  errors. Invalid RH refuses admission. This is an estimated ideal mixture:
  viscosity and conductivity retain the dry-air fits, the saturation
  enhancement factor is omitted, and no exact humid-medium prediction or
  experimental validation is claimed. G0/G3 assembly tests check the dry
  boundary, refusals, independent mixture density/speed/impedance, deterministic
  pressure, and waveform-derived duct pitch. The fixed eight-cell Cauer duct
  retains its first-order spatial boundary error; the pitch check uses its
  analytic LC frequency scale and does not establish 5% continuum accuracy.
  ISO absorption remains the sole
  bulk absorption inside its window; no extra classical term is added.
- String path, `EA = 0`: triangular-pluck (or bow) modal ICs marched by
  `ModalAcousticTimeModel`. The shared fixed-end string observer reads actual
  acceleration: odd modes give `p = rho A_odd q_ddot / (4 pi r)` in physical
  displacement coordinates. Even modes have zero monopole area and signed
  moment `M_even = -w L²/(n pi)`, giving `rho M_even q_dddot/(4 pi r c)`.
  This follows the spatial expansion of distributed retarded monopoles
  ([Wood, eq. 6](https://euphonics.org/11-7-1-compact-sound-sources-monopoles-dipoles-quadrupoles/)).
  Linear acceleration includes held forcing and viscous loss; nonlinear
  acceleration uses the admitted pHS storage and `(J-R) grad H + G u`, so
  stretching, conservative contact, and applied forces reach emitted pressure.
  Jerk is a backward difference of accepted accelerations: first-order at the
  reported endpoint, second-order at the interval midpoint for smooth motion.
  The initial acceleration is evaluated from the initial state and first held
  force; no zero-history startup spike or ideal release impulse is synthesized.
  The observer retains the linear path's 1e7 Pa pressure ceiling.
  This is an uncoupled compact strip approximation with the even-mode observer
  on the positive string axis. Common propagation delay, angular listener
  geometry, cylindrical-wire radiation, resolved impact impulses, and radiation
  back-reaction are not supplied by this observer. G1 tests check signed SI
  pressure against a damped oscillator, even-mode jerk and refinement, nonlinear
  stretching-force pressure, and silence under a held static load. Fletcher
  inharmonicity `ω_n = n ω_1 √(1+B n²)` when `EI > 0`. Stokes air
  drag from `GasState` plus the authored constant modal ratio. The reusable
  `air_path::oscillating_cylinder_air_resistance_per_length` computes the
  Desvages (2018) section 3.2.2 cylinder resistance approximation
  `R=2 pi mu (1+r sqrt(2 rho omega/mu))` [N s/m²]; both string paths use
  dimensionless `zeta=R/(2 mu_linear omega)`. Width is interpreted as circular
  diameter. The approximation retains a constant low-frequency continuation;
  added mass, end effects, confinement, rarefaction, turbulence and nonlinear
  fluid drag remain unclaimed. Explicit Rayleigh damping still replaces the
  complete loss sum when no material bending law is selected. The explicit
  Kelvin–Voigt path above instead combines only bending viscosity and air drag.
  Without a selected constitutive law, `damping_ratio` is a finite nonnegative
  constant per-mode reduced-model input. Zero explicitly assumes no internal
  loss; elastic EI does not create dissipation. No hidden `2e-7*omega` addend
  or relaxation spectrum inferred from a single ratio remains. The authored
  ratio is not a material measurement; physical loss requires admitted data
  and a selected Kelvin–Voigt or causal Prony law. G1 checks the third mode's
  actual SI pressure against the independently damped oscillator, and G3 checks
  equivalence to explicitly zero bending viscosity across linear, KC and Dirac
  paths. Invalid authored ratios refuse instead of being clamped.
  G1/G3 tests check the resistance formula and
  pressure ring-down under changed gas pressure. A second
  polarization at `1+detune` is a second fixed-end member on the same clock
  and shares every plate (it is not an independent unused body). Moving-end
  second polarization requires explicit transverse attachment geometry that
  this assembly does not carry; it refuses before realization, with or without
  Prony bending, instead of silently rendering only one plane.
- String path, `EA > 0`: `fs-nlmodal::kirchhoff_carrier_string`
  (fixed-fixed) or `kirchhoff_carrier_moving_end` (Dirac join)
  + Gonzalez
  pHS step, with the same stiff-string frequencies overwriting the
  linear ω_n. Bridge force `T_eff y'(0)` drives every listed body mode
  (Carcagno-band Sitka pair is a constructor, not a named guitar).
- Bow: a set-valued Stribeck contact (`stribeck_friction::StribeckLoadLine`).
  Each substep solves the exact-ZOH endpoint load line `v = v_free - m F`
  against the unregularized kinetic law of `fs-tribo::FrictionLaw::Stribeck`
  and its static cone `mu_s N`. The residual has at most one local maximum and
  one local minimum, so its stable roots are bracketed exactly. Where a stick
  reaction and sliding roots coexist, the contact keeps its branch
  (McIntyre–Schumacher–Woodhouse hysteresis) and changes only at a fold, a
  slip reversal or a lift. No capture window or ramp scale enters the bow
  path. G1 checks root selection against a dense-scan oracle, the hysteresis
  loop against an independent fold scan, and the Coulomb closed form. The
  bowed gates measure a 4..64 substep ladder at the playable point (playable
  on every rung, steady pitch within 12 cents). The former capture window
  (`max(stiction, 0.02)` and a 0.15 m/s flip basin) is gone; at 3.9 N, its old
  gate point, the regime was not converged in substeps under either contact.
  Helmholtz motion is possible with enough modes; it is not guaranteed and
  not a measured rosin curve. For explicit ODE embeddings elsewhere,
  `StribeckFriction` delegates traction to
  `fs-tribo::FrictionLaw::regularized_traction_1d`, reversing the port's
  driver-minus-body velocity to the owner's body-minus-driver convention. The
  shared owner supplies deterministic exponential evaluation and refuses
  malformed coefficients, negative/non-finite normal loads, and
  unrepresentable force. G1/G3 checks compare the force with an independent
  Stribeck formula and the work of equal-and-opposite contact forces. This
  wiring admits no sourced interface data or authored modal losses.
  `BowedStringCard::from_material_specimen` lowers a resolved circular specimen
  using the shared pinned-beam dispersion, Kelvin–Voigt loss and cylinder-air
  drag owners. It requires a matching loss/material identity and frequency
  band, rejects memory, moving supports and polarization splitting, and applies
  runtime bandwidth admission. The caller retains the source receipts; this
  fixed-state, small-amplitude card grants no new material authority.
  Elastic and viscous support shear are evaluated at the bridge, `x = 0`:
  `T y'(0) - EI y'''(0) - eta I y'''_t(0)`. Explicit cards must include the
  viscous bending contribution in their modal damping ratios. Bow velocity
  and bridge-force histories observe the same completed sample endpoint.
  G1 tests compare independent held-load motion, energy and support shear;
  G3 checks shared pluck/bow losses and sensitivity to material viscosity and
  ambient pressure. This does not evolve material state, route damping heat,
  model nonlinear stretch, or establish bow-contact passivity; substep
  convergence is measured at the playable point only.
  `FrictionIsland::InterfaceStribeck` consumes an ordered dry-interface card
  through `fs-tribo`, retaining its identity, history and caller authority.
  The complete sliding law, including velocity strengthening, drives force;
  other friction families refuse. Fixed temperature, load/nominal-area pressure
  and actual slip must remain inside the card domain at admission and both
  endpoints of each substep, including sticking. These sampled checks do not
  bound interior extrema. G1/G3 compare signed, load-dependent motion to an
  independent oscillator and check runtime domain departure. This path supplies
  neither measured-data admission nor heating, wear or evolving contact area.
  Bowed plate outputs report signed volume velocity in m³/s, not point velocity.
  Plate observers admit positive finite listener distance before stepping and
  refuse non-finite volume flow or pressure. G0/G3 checks cover malformed
  distances and inverse-distance pressure scaling with identical motion.
  This compact observer adds no propagation delay or general near-field claim.
  The material string bank can retain both transverse displacement planes
  through `MaterialStringSource::with_transverse_polarizations(TransverseDirections::try_new(bow, bridge)?)`.
  Both planes use the same original circular section, installed tension, density,
  bending stiffness, material loss and gas drag. Their physical directions project
  displacement and force with the same work-conjugate map; shared body compliance
  determines polarization mixing without an authored detuning. Coordinate order is
  first transverse plane, second plane, then optional material torsion. Surface
  bow/pluck ports include the original radius-times-angle contribution; centreline
  loads and unilateral stops carry no invented radial torque. Stops follow the
  source's declared centreline/bow axis. Existing explicit torque/rotation bridge
  connections compose with both planes and retain their own source-derived angular
  maps and passive laws. Material body/cavity pressure and scheduling remain in the
  original joint solver.

  `frankensim-material-strings-v3` extends the v2 section/support grammar with one
  required `polarizations single` or `polarizations two BOW_X BOW_Y BRIDGE_X BRIDGE_Y`
  record after each string section. Both directions must be finite unit vectors.
  Existing v1/v2 files and optional `rotational-bridges` records retain their
  semantics. `polarized-pluck.material-strings` exercises the data-to-PCM route.
  This circular, frozen-state linear model does not supply arbitrary
  six-degree-of-freedom bridge motion, independent stop-normal geometry, shared
  nonlinear two-plane tension, helical winding/slip, or a measured-instrument timbre
  guarantee. Nonlinear material source conversion uses the existing shared
  two-plane Kirchhoff–Carrier stretching law.

  Seven focused physical tests cover material-mode identity, rotational
  covariance, directional sound changes, force work conjugacy, constructor
  ordering/refusals, combined rotational-bridge response and two-plane unilateral
  contact. Three example tests cover source input, material/angle swaps,
  PCM/callback identity and admission. Native verification is recorded separately.
  `ReciprocalStringBank::with_cylinder_gas_loading` selects a causal gas
  operator for the same circular material strings. The shared
  `air_path::cylinder` owner derives the unsteady incompressible Stokes impedance
  `Z(s)=pi*rho*a^2*s+4*pi*mu*z*K1(z)/K0(z)`, `z=a*sqrt(rho*s/mu)`, from
  actual outer radius and the retained complete gas. The continuum follows
  [NASA/CR-2007-215050, equation (7)](https://ntrs.nasa.gov/api/citations/20080006726/downloads/20080006726.pdf).
  Its positive spectrum per log relaxation rate,
  `8*mu/[pi*(J0^2+Y0^2)]`, supplies physical free inertias joined by passive
  dashpots in the original exact reciprocal system. Both transverse planes
  carry instantaneous added gas mass and viscous memory; original solid
  density, bending, installed tension, material viscosity and torsional inertia
  remain distinct. The active transverse card replaces sampled gas resistance
  with intrinsic material loss plus the explicitly collapsed slow tail, so
  fluid damping is not counted twice. The same intrinsic-loss owner supplies
  Kelvin–Voigt bending and the construction's independent tensile loss factor
  in the loaded mass basis; selecting causal fluid loading must not erase
  either. Prony bending keeps its internal strains alongside that tensile
  loss, without inventing a second bending-viscosity law.

  All material-string performance versions accept the optional
  `cylinder-gas-loading X_MIN X_MAX INTERVALS` row before an optional
  `prepared-linear` row and then `cavity`, after any section, polarization and
  rotational/translational-bridge records. The same
  original geometry/gas supplies the operator; no extra material coefficients
  are authored. This binding precedes additional bows, stop profiles, static
  preload and the physical pressure/PCM schedule. File construction counts all
  mechanical/fluid coordinates and string/body/torque/cavity/tensor-loss ports
  against the explicit mode, connection and complete setup-work allowances.
  The source fixtures and regressions compose material/mixture substitutions,
  regional body loss, both transverse planes, torsion, cavity, stops and release.

  These fluid components append after physical bodies. Full tensor body loss,
  cavity connections, pressure observation, dry interface and torque maps
  survive binding, either cold torque/fluid order, and multiple-bow scheduling.
  Pressure and force maps follow the changed transverse mass normalization;
  fluid coordinates receive no fictitious direct acoustic pressure weights.
  Added kinetic storage, reciprocal reactions and viscous heat enter the same
  complete work balance. Static initialization admits exactly zero-stiffness
  positions only when their applied static force is exactly zero. Their
  already-zero positions remain gauges while their velocities remain physical
  free inertias; no artificial spring or general singular preload solve enters.

  `CylinderFlowResolution` declares independent slow/fast dimensionless
  wavenumber cutoffs and log-grid intervals: `0<X_MIN<=0.5`,
  `X_MIN<X_MAX<=3000`, and `1..=4096` intervals. The exact integrated slow
  tail is collapsed to resistance with dynamic error bounded by
  `R_low*lambda_min/omega` at positive frequency. The final log cell and the
  complete faster continuum share one positive mass/dashpot branch. For that
  original spectral measure, `M=integral(dR/lambda)` and
  `N=integral(dR/lambda^2)` give `lambda_eff=M/N`, `R_eff=M^2/N`. This retains
  the added inertia and leading dissipative response: both exact and reduced
  impedances start with `s*M-s^2*N`, and their difference starts at order
  `s^3`. The fast-tail reduction applies to transverse, axial and rotary flow;
  instantaneous potential-flow mass and the rotary steady resistance remain
  separate. Fluid inertia is evaluated from density and geometry, so a
  viscosity-only change changes resistance/rates without changing storage.

  This completion keeps exactly `INTERVALS` memory coordinates, with no extra
  runtime operations or new material parameters. It adds cold positive-weight
  moment quadrature: logarithmic panels below the Bessel crossover and inverse
  wavenumber above it, with cancellation between refinements. The large-
  argument squared Bessel modulus uses phase-free Hankel amplitudes rather
  than extrapolating the individual oscillatory J/Y phase certificate. The
  internal moment-quadrature tolerance is `1e-10`; it does not certify a
  finite-frequency impedance or coupled trajectory. In particular a finite
  number of branches cannot reproduce the continuum's unbounded square-root
  high-frequency law. The upper cutoff must still be refined for the physical
  bandwidth; the complete fast tail is represented by its first two inverse-
  rate moments, not resolved at every rate.
  Refine all three axes and the actual coupled response. Every branch adds one
  coordinate and one connection per retained transverse/twist mode, and per
  axial mode when longitudinal waves are selected. The shared mechanical
  connection ceiling is 4096, with caller mode/connection/setup budgets still
  enforced before allocation; that port capacity alone does not raise the
  ordinary exact connected-island limit. An explicit prepared linear budget
  admits the complete larger system as described below. Repeated/live bindings and loaded elastic
  modes outside their original material bands refuse.

  Operator regressions compare complex force with an independent modified
  Bessel reference and exercise every cutoff/grid axis. Material-bank tests
  exercise pressure-dependent inertia, viscosity-only changes in actual body
  pressure, both transverse planes, unchanged torsion, fluid kinetic energy,
  torque binding order and preload/release. Source review is not native test
  execution. This local long-cylinder, small-motion Newtonian incompressible
  model supplies no finite-amplitude convection/wakes, axial-gradient/end/confinement
  or neighbouring-string interaction, rarefaction, compressible
  radiation, evolving atmosphere or measured violin/cello/bass fidelity claim.

- Reed path: quasistatic or massive Bernoulli valve. Reeds with explicit
  positive stiffness derive their effective pressure area from the
  declared closing relation `A = k H / Pc`. Characteristic-line motion and
  ODE pressure force/swept-face flow share that area; a fixed 25 mm face
  length no longer overrides explicit stiffness and closing pressure.
  This is an effective modal area, not measured geometry. The zero-stiffness
  legacy reduction still derives stiffness using a 25 mm face length, and
  damping is now an explicit finite nonnegative ratio, consumed as
  `c = 2 zeta sqrt(k m)`. Reduced valve cards pass `zeta = eta/2`, preserving
  their existing `c = eta sqrt(k m)` model. This removes the runtime's fixed
  0.35 substitution; authored fixture ratios remain explicit assumptions,
  not measurements. G1 covers actual velocity decay; G4 rejects invalid ratios.
  G1 tests exercise the force balance and actual moving-reed acceleration.
  The massive characteristic-line junction now uses the current opening and
  swept-face flow `-A v_mid`, plus body volume flow, in the pressure solve. Positive
  velocity widens the opening: the ODE and characteristic paths share this
  sign, conjugate to `F = -A (Pmouth - Pbore)`. It corrects the ODE's former
  positive swept-flow sign; see Silva et al., arXiv:0810.2870, Eq. 27.
  Its fixed-opening residual is monotone, bracketed between zero jet flow and
  zero pressure drop, and uses the same signed Bernoulli primitive as the ODE.
  The massless solver and its Newton Jacobian are not used for massive reeds;
  fast-mode massive samples count as strict fallbacks. G1 tests cover closed
  apertures, reverse flow, and instantaneous pressure/face work cancellation;
  G3 covers mass-dependent bore pressure at equal static compliance.
  Characteristic-line linear mechanics now uses midpoint position and velocity,
  with that same average velocity in swept flow. Eliminating the linear reed
  response produces a positive effective load for the existing aperture solve.
  G1 checks the discrete spring/kinetic energy change against pressure work,
  Bernoulli loss and viscous loss without contact, plus energy preservation and
  second-order phase convergence of the isolated oscillator.
  Aperture geometry and elastic lay force are held at the start of the step.
  Single-slit Hunt–Crossley loss uses the same midpoint velocity as mechanics:
  the active reaction is `elastic - elastic*chi*v_mid`; a second zero-reaction
  solve handles the nonadhesive unloading clamp. The coefficients come from
  the admitted fs-dcontact obstacle, restricted to one unit opening coordinate.
  G1 exercises loading, unloading and within-step reversal; a lagged-velocity
  negative control injects energy where the accepted loss removes it.
  This proves the sign of loss work against the held elastic reaction, not
  that the latter matches the change in nonlinear contact potential.
  The isolated ODE reed retains its own pHS stepping scheme. This is not a
  nonlinear contact or whole-bore energy certificate: contact work, coupled
  geometry/time refinement and onset validation remain MR68 obligations.
  Closure stop (MB 3381 and 3605, 2026-11-05). The reed coordinate cannot go below zero opening: the
  reed cannot pass through its mouthpiece. After the implicit midpoint solve above, which is unchanged,
  a step that starts at `y >= 0` and whose solution ends at `y < 0` is replaced by the unilateral
  constraint `0 <= y` with a non-negative multiplier: the reed ends the step at exactly 0 with mean
  velocity `-y/dt` and end velocity 0 (rigid and inelastic), and the aperture junction and the lay
  contact are evaluated at that landing velocity by the same code. The multiplier is the stop's push,
  the momentum residual at the landing velocity; the residual rises with the mean velocity, so it is not
  negative, and a negative one refuses ("closure stop would pull the reed shut"). Flow through the shut
  channel is the existing `max(0, midpoint opening)` clamp: zero while the reed rests shut; only the
  landing step, in which the opening is still closing, passes the flow of its half-open midpoint.
  A step whose free solution stays at or above zero is the original arithmetic, bit for bit
  (`step_with_contact_response_and_stop` with `closure_stop = false` is that original step). A reed that
  enters a step already below zero is left to the free solve. The loaded finite-gas path
  (`reed_bore::loaded`) takes no stop: its work ledger balances the heats and works of the free step and
  has no impact term, so a stopped step would refuse there. The plastic stop is a modelling choice: no
  restitution is printed for a reed on its lay. Measured on the clarinet G4 test card without its slit
  line (impact about 0.5 m/s, about 5e-7 J per impact, 0.2 to 0.4 mW): at 3.0 kPa the earlier law left
  the opening below zero on 14.75 % of the steady samples (minimum -0.0226 mm) and at 3.6 kPa on 38.22 %
  (-0.0450 mm); with the stop the minimum is 0 and the reed is shut for 20 % and 36 % of the samples;
  H8 to H12 over the recording moved by -0.4 and -2.5 dB. The live clarinet cards at 1.8 kPa never close
  (minimum opening 0.080 mm at G4) and are byte-identical (G4, D3, Eb5). Other families move, as they should: RubyPrairie's
  alto-sax cards were below zero on 18.0 % (C4), 17.1 % (E3) and 19.0 % (F3) of the steady samples (minimum -0.96,
  -1.25 and -1.29 mm), a held C4 card with the clarinet lay on 8.0 % (-0.18 mm), and the oboe Bb4 noise-row card on
  19.6 % (-0.30 mm); those renders change and are re-rendered on the landed binary. No bassoon card was held.
  G1 tests: the free reed pressed through a soft lay
  stops at exactly zero opening and stays shut with zero velocity, and is bitwise the free step until
  the crossing; a reed that stays open is bitwise unchanged; the flow of the shut channel is zero and
  the junction wave is the closed one; a stop that would pull the reed shut refuses.
  Isolated
  cylindrical bores use the `acoustic_chain` ODE with a
  `ViscothermalPin` (massive
  reed = `mass_spring_damper` + jet + face flow). Blow/reed on a
  moving-end string×plate×duct is a leftover `join_port` inlet.
  Moving-end string×duct without a plate is
  `transformer(waveguide, chain, A_inlet)`. Fixed-fixed
  string×duct without a plate shares the ODE clock without
  a force port (`φ(0)=0`).
  Reed lay
  is `fs-dcontact` (not a private Hunt–Crossley).
- `ThinPlate` modes come from `fs-plate` + `fs-modal` (certified
  eigenpairs), not caller Hertz. Mouth pressure × area drives the
  plate; plate volume velocity returns to the waveguide
  (structure–bore).
- Fixed-end Kirchhoff–Carrier members prepare the complete original
  contact-wrapped, bending-memory system at the render clock.
  `PreparedNonlinearZoh` advances its quadratic part exactly and solves the
  original nonlinear discrete-gradient endpoint. Admission uses the actual
  `2 * state_dim + port_dim + 1` propagation and `2 * state_dim` heat matrices
  against the unchanged 256-coordinate ceiling. The held nonlinear-force
  shadow reuses the original admitted Q/J/R spectral data. For bare KC this
  permits 51 complete modes; attached memory consumes the same state budget.
  The native account retains original Hamiltonian change, independently
  integrated physical resistance heat and same-window held-force work.
  Candidate state, force and mean-output buffers are reused. Every model
  admitted by the former midpoint size reservation fits this prepared path,
  so the redundant private fallback is removed. Complete models above the
  actual capacity refuse, and preparation failures propagate.
  Every selected mode, conservative obstacle and memory arm remains present.
  Retained bare/bending frequencies must also fit the linear modal observer's
  Nyquist guard. This guard does not certify nonlinear or contact harmonics:
  time, modal and contact refinement remain separate accuracy obligations.
  A changed clock refuses before advancing any history. Failed KC trials
  restore the driven texture path and retain accepted mechanics/observer
  history; the surrounding multi-member render is not a joint transaction.
  Bridge-force outputs and compact pressure use the same accepted coordinates
  and actual force/acceleration maps. Fixed-end support coupling remains the
  original one-way force path; moving-end descriptor stepping is separate.
  Focused `acoustic_realize::kc_step_tests` exercise the high-frequency
  quadratic limit, nonlinear period refinement, complete contact/memory work,
  sampled-band refusal and rejected trials. Added regressions retain the
  highest partial of a 43-mode bank beyond the former split capacity and
  refuse a complete 52-mode bank without truncation. Native execution of
  these additional regressions remains pending.
- String obstacles are `fs-dcontact` power-law potentials on sine
  collocation. On the Kirchhoff–Carrier path the conservative
  potential is `ContactStorage` wrapping the KC Hamiltonian, so
  the prepared step retains its contact energy. Optional `mu_kinetic` adds
  `fs-tribo` Coulomb traction as a modal port force (not a gradient
  of `H`). The linear modal path still applies both as a force.
  Hunt–Crossley `χ` is a scenario field on the obstacle and a
  dissipative port force (`modal_hunt_crossley_forces`), never a
  term in `H`.
- Reed and structure–bore time port: TMM `R(ω)` sampled on the DFT
  grid and inverse-transformed to an FIR. Isolated linear blow (no
  body) is the same `DelayedFilter` object filled with `IFFT[Z_in]`
  (or the mouth transfer `p/u` when the end is open). A vented
  reflectance FIR does not ring a measurable period; the impedance
  FIR does, so tone-hole shortening is TMM-emergent on one port
  type. There is no one-pole `TravelingWaveLine` fallback on those
  paths. Unflanged-open Nyquist `ka > 1` is refused; a
  flanged mouth uses the Rayleigh piston and is not.
- Card-driven acoustic compilation: `material_assembly::compile_material_assembly`
  resolves exact `MaterialCard`/`QueryPoint` inputs with explicit claim selection,
  then calls the existing circular-string, uniform-plate, regional-chart and
  local duct-wall binders. Only consumer-needed scalar properties are resolved: rho/E
  for elastic strings, the selected in-plane plate constants, and additional
  viscosity or equilibrium/modulus/time fields for explicitly selected sourced
  string/plate bending loss and explicitly selected isotropic plate thermal loss.
  Duct walls use the same source selectors as graph and DG wind performances,
  including isotropic or principal-aligned elasticity and selected material memory.
  Inactive property families are not queried. Loss selectors
  are retained; old numerical loss coefficients and citations are rebuilt.
  `StringPrestressBinding::ThermalExtension` selects the existing alpha(T)
  integrator and sourced thermal string binder. Its reference point and alpha
  claim selection are explicit; the immutable card and current point come from
  the same string source as density, elasticity and selected loss properties.
  Mechanical/loss pins and expansion pins cover their separate requirements.
  Piecewise-linear expansion curves are integrated across all crossed knots;
  current alpha times temperature change is not substituted for that integral.
  Fixed supports use the reference stress-free length to derive elastic strain
  and tension. Fixed mass or current radius remains an independent constraint.
  The specimen retains both mechanical and integrated-expansion receipts, and
  the ordinary acoustic realizer uses the resulting tension and rebuilt loss.
  `Prescribed` keeps explicit force/extension/tuning constraints and raw authored
  eigenstrain; it does not query or claim sourced expansion. Invalid paths,
  missing/ambiguous data and slack/compressive states refuse through the existing
  owners. No thermal transport, finite deformation or melting is added.
  `CompiledMaterialAssembly` keeps the numeric assembly and source-bound specimens
  immutable, and realization calls the existing ordinary/chart dispatcher.
  A source-bound duct wall requires explicit physical receiver geometry through
  `realize_with_options`; it cannot silently disappear into an unsupported path.
  Unbound components remain authored. Material conditions remain separate from
  ambient-gas conditions; this compiler does not infer thermal equilibrium or
  upgrade the authority of the whole assembly. Geometry, constraints, forcing,
  supports, damping choices and solver controls retain their existing meanings.
  Resolution/binding errors preserve component context and the owning diagnostic;
  failure leaves the input and prior compiled specimens intact. A thermal plate
  selects the existing six-property thermoelastic resolver; optional bending
  viscosity is requested in that same atomic bundle, including exact claim pins.
  The resulting state rebuilds the fixed-mass/thickness specimen and both loss
  mechanisms. Stale numerical coefficients and citations are replaced together.
  Regional bindings explicitly select `RegionalPlateLoss::Authored` or
  `KelvinVoigt`. The latter resolves every region's loss with its elastic state:
  isotropic regions require uniaxial bending viscosity, while orthotropic regions
  require all six reduced plane-stress viscosity coefficients, including signed
  or explicitly zero couplings. The original tensor owner handles rotation and
  admissibility. `realize_assembly_with_material_plate_chart_and_options` keeps
  the resolved chart through the common acoustic dispatcher, so its existing
  modal curvature-energy loss projection and frequency-domain refusals cannot
  disappear at a bare-chart conversion. Authored damping must be zero alongside
  sourced loss. This retains the existing diagonal weak-loss image, not full
  nonproportional damping. G1/G3 compare uniform eta/E against its analytic modal
  decay, substitute one region's loss source, and observe changed coupled
  pressure and free decay without changing elastic frequencies or the gas state.
  An orthotropic material selection with this isotropic thermal law refuses.
  Specimen temperature must be an explicit positive absolute T coordinate; the
  surrounding gas does not supply it implicitly. Compilation is
  synchronous with no cancellation guarantee; ordinary runtime admission runs at
  realization. This is a Rust description API, without a new wire format, global
  scenario identity, optical binding or conservative material evolution.
  G0/G1/G3 exercise explicit-binding pressure parity, consumer-minimal exact-source
  retention, simultaneous string/plate card substitution at fixed mass, regional
  chart substitution through moving-end string/cavity pressure, active Kelvin–Voigt
  and Prony coefficient rebuilding, and missing-member/property/domain/law refusals.
  Thermal compiler tests vary specimen temperature through source curves and
  check changed pressure in linear, sine von Karman and sampled-FE plate paths,
  with and without separately selected viscosity. A fixed-mass card swap checks
  rebuilt thickness and both damping addends against independent slab diffusion,
  mean-curvature participation and Kelvin–Voigt arithmetic; conflicting
  heat-capacity claims require the selected
  pin to survive compilation. These fixtures validate implementation, not an
  experiment or a transient heating simulation.
  String thermal compiler tests independently calculate integrated strain,
  fixed-mass area/inertia, tension and stiff-string frequency, then measure
  emitted-pressure pitch at two temperatures with and without sourced viscosity.
  The viscosity fixtures supply scalar coefficients at explicitly supported
  temperature states, preserving the existing frequency-wide constant-law
  requirement; sampled coefficients are not promoted to broadband loss laws.
  They also check exact alpha selection independent of mechanical pins, force
  control without alpha resolution, and unchanged inputs after path/slack refusals.
- Uniform plate material binding: `with_uniform_plate_material_state` consumes
  the same `ResolvedMaterialStatePoint` carrier as the circular-string and
  disc adapters. Isotropic elasticity requires rho/E/nu; principal-axis
  orthotropic plane stress requires rho/E1/E2/nu12/G12. It retains every original
  property receipt and checks exact scalar quantity schemas. Card axes map to
  section material axes, with reciprocal nu21 on exchange. The explicit
  `ThinPlate::material_angle_rad` rotates that uniform section counterclockwise
  from rectangle x, including D16/D26 bending/twisting coupling. Geometry,
  supports and forcing remain fixed. The angle participates in the v2 specimen
  identity. No out-of-plane constants are invented in either
  linear or FE-sampled acoustic plate construction. A fixed-thickness or
  fixed-mass prescription derives thickness and mass coherently; geometry and
  material/model identity are retained separately from supports, pretension,
  authored damping and solver controls. These are independent frozen specimens,
  not a conservative thermal evolution. The purely elastic binder refuses stale
  thermoelastic addends. `with_uniform_isotropic_thermoelastic_material_state`
  instead rebuilds the mechanics and loss inputs from the existing resolved
  thermal bundle, retaining its complete receipts and optional viscosity data.
  The source-bound orthotropic adapter also refuses geometric
  nonlinearity: the existing sampled nonlinear plate's isotropic membrane
  approximation is not a resolved orthotropic membrane law. Anisotropic thermal
  loss and phase changes remain unavailable here. The raw
  nonlinear path checks E, nu and G before selecting isotropic sine modes.
  Its rotated FE bending still uses the documented isotropic membrane
  approximation. The modal search uses a sine trial Rayleigh quotient;
  it is not an exact off-axis frequency or a completeness certificate.
  G0/G1/G3 cover receipts, reciprocity, fixed-mass geometry,
  semantic aliases, model boundaries and real pressure changes from material
  and axis swaps, with aligned eigenfrequencies checked against Navier theory.
  An arbitrary angle changes actual assembly pressure; a quarter-turn recovers
  the pressure and eigenfrequency of reciprocal axis exchange.
- Uniform material-bound plate viscosity: `ResolvedPlateSpecimen::with_kelvin_voigt_bending_loss`
  selects proportional isotropic Kelvin–Voigt bending. The dimension-only
  `kelvin_voigt_bending_viscosity` [Pa s] claim must be nonnegative and carry an
  explicit frequency band. Rho, E, nu and eta must be validity-wide scalar
  constants or temperature-only curves at fixed positive absolute T; their
  frequency bands intersect. Temperature, moisture and other
  coordinates remain at the admitted state point, with original uncertainty
  and observation status retained. The assumed viscous plane-stress tensor has
  the same Poisson ratio as the elastic tensor: the moment law is
  `M = D (curvature + (eta/E) curvature_rate)`, giving `C_b = (eta/E) K_b`.
  This is a proportional constitutive specialization of
  [Kelvin–Voigt strain-rate damping](https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.029.html#3441080),
  not a measured general bulk/shear viscosity tensor.
  The ordinary linear, sine von Karman and sampled-FE plate paths consume it.
  Each retained reference-mode frequency must lie in the admitted band, even
  for a zero-strength viscosity claim. Only bending energy receives this loss;
  geometric prestress does not. For FE modes the diagonal coefficient is
  `(eta/E) phi^T K_b phi / (phi^T M phi)`. Off-diagonal damping coupling is
  omitted, a modal approximation when prestress makes damping nonproportional.
  The nonlinear path retains bending-rate damping but supplies no nonlinear
  membrane viscosity or guarantee that finite-amplitude frequency shifts and
  generated harmonics remain in band. Authored scalar damping alongside this
  law refuses. Explicit thermal loss remains separately additive; its inputs
  must exclude this viscosity's contribution to avoid double counting.
  `compile_material_assembly` resolves the selected viscosity with the plate's
  elasticity and density. A material rebind recomputes its coefficients and
  source identity rather than carrying the old material's viscosity forward;
  missing source data refuses. Thermal rebinding of an already viscous raw
  plate refuses because it would otherwise retain stale viscosity inputs.
  G1 checks the actual sine-mode damping operator and ringdown against the
  analytic Kelvin–Voigt plate equation with timestep refinement. Independent
  geometric-energy integration checks prestressed DKT and sampled-FE damping,
  including additive thermal loss. G3 changes only viscosity through the card
  compiler, retaining elastic eigenfrequencies while changing actual assembly
  pressure in all three plate paths. Fixed-temperature curve tests retain the
  source interpolation receipts and compare plate frequencies/damping, every
  Prony string branch and actual pressure against independently evaluated
  coefficients. The thermal string compiler also exercises interior T values
  with simultaneous integrated expansion and bending viscosity. These are
  frozen-state calculations, not constitutive evolution during heating.
  These synthetic implementation fixtures
  are not experimental material/ringdown validation.
- Regional linear plates: `ResolvedPlateSpecimen::section` exposes the rotated
  numerical section for assignment to `PlateChart` elements. The caller retains
  the resolved specimens for source receipts; a numerical chart does not mint
  material authority or a regional provenance identity. The full-rectangle
  specimen mass is not reused as regional mass: fs-plate integrates local rho*h
  and rotary inertia from each actual section.
  `compile_plate_material_chart` additionally binds disjoint named regions to
  exact resolved material states, constitutive/axis choices and fixed-thickness
  or fixed-region-mass prescriptions. The immutable `ResolvedPlateChart` retains
  these inputs alongside its numerical chart. It integrates actual triangle
  areas to derive h = m/(rho A_region); `with_region_material` rebuilds sections,
  mass and identity using the retained constraints. Failure leaves the original
  binding unchanged. Empty/duplicate names, overlapping or missing assignments
  and out-of-range triangles return `PlateAssignment`; geometry and section
  refusals preserve the owning `PlateError`. No fallback material is inserted.
  The specimen identity includes mesh geometry/connectivity, boundary-node set
  and per-element source/model/thickness/angle in mesh order. Names, assignment
  order and repartitioning with identical resolved per-element inputs do not
  alter it. It is not a complete acoustic-scenario or authored-constraint identity.
  Compilation is synchronous, with no cancellation guarantee or material-state
  evolution; cloned numerical charts do not automatically retain source authority.
  `certified_chart_radiators` consumes that chart with `PlateChartRadiation`
  (explicit eigenvalue search interval, mode budget, support/prestress, authored
  damping and nodal unit-force footprint). All force weights must be finite
  and sum to one; supported forces react into the support without redistribution.
  `realize_assembly_with_plate_chart` feeds the resulting modes into the same
  assembly dispatcher, body ports, piston loads, string/cavity/duct coupling
  and pressure observer as uniform plates. A simultaneous `assembly.plate`
  refuses. Extra authored lumped bodies still compose. The modal interval
  need not include the fundamental; certification applies only in that window.
  It does not certify stability or exclude buckling modes outside that window.
  `ResolvedPlateChart::with_kelvin_voigt_bending_loss` selects regional scalar
  loss (`RegionalBendingViscosity`): isotropic regions use the proportional law
  above; orthotropic regions use their six explicit `plane_stress_viscosity_*`
  coefficients, resolved and rotated by the same owner as the full tensor image
  (`plate_viscosity::resolve_region_viscosity`; no coefficient is inferred). The
  one diagonal evaluator, `plate_viscosity::DiagonalPlateViscosity`, gives
  `zeta_i = sum_e w_r q_e^T K_e q_e / phi^T K phi`: isotropic regions use the
  elastic DKT matrix with `w = tau omega / 2` (the legacy arithmetic, moved bit
  for bit); orthotropic regions use the DKT matrix of the rotated viscous tensor
  with `w = omega / 2`, equal to `C_ii / (2 omega_i)` of the full projection for
  an exact eigenpair. This image omits off-diagonal modal damping;
  `PlateViscosity` projects the full tensor, and a chart with scalar loss
  selected cannot also select it. Tests:
  `isotropic_regional_loss_routed_through_plate_viscosity_is_bitwise_legacy`
  (zeta and radiator bits against a verbatim pre-move copy) and
  `orthotropic_region_now_admitted_and_matches_full_matrix_diagonal` (1e-13 plus
  twice the measured eigenpair residual).
  No spatial thermal-loss law, nonlinear membrane, remeshing, phase evolution
  or full-wave radiation is claimed. Interfaces are perfectly bonded, with
  thickness centered on one flat mid-surface as specified by fs-plate.
  G3 checks uniform/element partition pressure equivalence, legacy uniform
  path parity including moving-end/cavity ports, and pressure/mode changes
  from local density, stiffness, thickness and grain on an unstructured chart.
  G0/G1/G3 additionally check assignment refusals, exact source retention,
  integrated regional mass, region-order/partition pressure equivalence and
  actual pressure changes after named material replacement under both mass
  and thickness constraints. These use synthetic materials, not experimentally
  validated material data or measured sound.
  G0 checks force footprints, modal controls and conflicting descriptions.
- Linear plate projection: the eigenvector, its modal mass `phi^T M phi`,
  signed radiating area `integral(phi dA)`, and signed bridge-force projection
  remain in one basis. Surface integrals use triangle-area P1 nodal quadrature;
  the declared drive footprint is uniform traction over the generated
  rectangle's first cell strip, normalized to one unit of total force.
  There is no mass inferred from drive amplitude, absolute-value radiation,
  or artificial area floor for cancelling modes. Generalized modal force is
  `b_drive F + A_mode p_face`. The pHS ports divide both coefficients by
  `sqrt(modal_mass)`, and its observers convert acceleration back by the same
  factor. Thus eigenvector sign/scale cannot change physical pressure or work.
  This follows modal virtual work and the compact limit of the
  [Rayleigh integral](https://euphonics.org/4-3-2-the-rayleigh-integral-and-the-baffled-piston/).
  The optional diagonal piston self-load uses the physical panel area to set
  its equivalent circular radius; input velocity is `A_mode q_dot/A_panel`
  and generalized reaction is `-A_mode p_face`. Modal normalization cannot
  change that geometry. The approximation omits mutual radiation impedances,
  higher multipoles and higher-order DKT trace reconstruction; eigenvalue
  certification does not certify radiation accuracy. G1/G3 check manufactured
  quadrature/mass, zero-monopole cancellation, mode sign/scale invariance with
  and without the actual fitted self-load, conjugate pHS ports, and density
  scaling through the actual assembly pressure output.
- Explicit anisotropic thermoelastic plate bending is available through
  `AnisotropicPlateThermal`, `certified_radiators_with_anisotropic_thermoelastic`,
  and `VkBody::from_plate_with_anisotropic_thermoelastic`. Thermal expansion is
  `[alpha11, alpha22, 2 alpha12]` in the original section's material axes;
  temperature, heat capacity at constant stress and through-thickness
  conductivity are supplied explicitly. The original oriented plane-stress
  tensor and DKT curvature retain the rank-one `Q alpha` thermal channel,
  signed cross terms, lossless curvature directions and prestress dilution
  through actual total modal strain energy. The material owner retains the
  complete insulated-slab diffusion spectrum without a thermal-mode loop.
  Both consumers keep the existing diagonal weak-loss modal reduction at
  nominal elastic frequencies, for a homogeneous frozen specimen state.
  Thermoelastic stiffness dispersion, strongly coupled modal clusters,
  lateral or boundary heat transport and transient thermal history are not
  solved. Caller-authored constants do not acquire material-pack authority.
  Duplicate thermoelastic loss and nonzero authored modal damping refuse;
  existing viscosity and material-source restrictions remain in force.
- Full thermoelastic memory is available through
  `VkBody::from_plate_with_thermoelastic_memory` and
  `thin_plate::thermoelastic_linear::ThermoelasticPlateImage`. Both consume the
  same explicit `AnisotropicPlateThermal` data and
  `ThermoelasticMemoryOptions`. The coupled material owner converts
  constant-stress to fixed-strain heat capacity; isothermal elastic stiffness
  remains the equilibrium stiffness. Original DKT curvature factors are
  projected through the same retained eigenvectors and coordinate scales as
  the mechanical body. Power-of-two scaled QR compresses spatial factors
  without a rank tolerance, retaining signed cross-mode coupling and weak
  independent channels. Positive diffusion branches carry their own storage,
  heat and histories through the existing pHS/held-force Maxwell owners.
  No parallel thermal damping ratio is inserted.
  The nonlinear sampled route retains the complete original requested basis
  under its existing pHS admission; `PlateBankImage` carries the complete
  appended histories. The linear image
  supports up to 64 requested modes, retains the original mass, force and
  reciprocal volume maps, fills drive buffers without allocation and prepares
  the existing `CoupledModalSystem`. Preparation keeps its original resource,
  energy and instantaneous-stiffness Nyquist admission.
  Explicit thermal-term, memory-state and projection-work allowances bound
  setup. Both consumers admit the represented instantaneous linear stiffness
  against the finite constitutive band; this does not bandlimit forcing or
  nonlinear response. The reported spectrum error covers only normalized
  complex constitutive response in that band, excluding spatial/modal
  truncation and projection/runtime arithmetic. Homogeneous frozen-state,
  insulated-face and through-thickness conduction assumptions remain;
  finite-amplitude membrane thermal loss, lateral/boundary heat flow and
  evolving mean temperature are not modeled. Existing duplicate-loss and
  source/viscosity refusals remain active; caller data gain no material-pack
  authority.
- Plate damping: a single authored viscous ratio is the two-point
  Rayleigh fit (`fs-material::visco`) through that ratio at ω₀ and
  4ω₀ so higher certified modes sit on the stiffness limb. Both linear and
  von Karman banks add thermoelastic bending loss only from explicit specimen
  inputs. `with_isotropic_thermoelastic_state` binds a resolved six-property
  bundle to the plate, moving elasticity, density, thermal properties, and
  temperature together and retaining its state identity. There is no
  density-threshold material guess or fixed-temperature substitute. Missing
  thermal inputs leave that loss unavailable/unclaimed; requesting the scalar
  law on an orthotropic plate refuses (including a mismatched shear modulus
  when E1=E2). The material owner evaluates the complete through-thickness
  thermal diffusion spectrum of a homogeneous isotropic slab with insulated
  faces. That scalar entry does not establish an anisotropic loss law or
  solve a transient thermal field. Callers must exclude this addend
  when authored damping or resolved thermal transport already includes it.
  Thermal participation is `E h^3/[12(1-nu)^2] integral (kappa_xx+kappa_yy)^2`
  divided by total modal stiffness. Original DKT quadrature gives a sum of
  squares, so saddle bending and twist create no spurious thermal loss.
  Sine modes use the plane-stress factor `(1+nu)/(1-nu)` and prestress
  dilution `D k^2 / (D k^2 + T)`. Conservative prestress enters only the
  total energy denominator; viscosity retains its own bending/total ratio.
  Authored loss is unchanged. The same projection reaches the sampled-FE
  von Karman bank. Curvature coupling and the diffusion spectrum follow
  [Norris and Photiadis, eqs. (33)–(45)](https://arxiv.org/abs/cond-mat/0405323).
  This is weak-loss, diagonal modal bending damping about the prestressed
  equilibrium; it does not model finite-amplitude membrane thermal loss or
  off-diagonal thermal coupling. The energy separation follows
  [Fedorov et al., Eq. (1) and Sec. II](https://arxiv.org/abs/1807.07086).
  G1 tests independently integrate geometric energy and compare actual
  unforced sine-mode ringdown to the analytic damped oscillator under refinement.
  G0/G1 tests inspect both assembled damping operators against independent
  thermal-mode sums and tensor-energy integrals, vary equal-density specimens
  and temperature, and cross
  the former density threshold without switching laws. Radiation
  reaction on linear compact radiators uses the sole Bessel–Struve piston
  owner in `fs-phs`, converting acoustic p/U to face p/v by the physical
  area, then fitting and passivity-repairing with `fs-vfit`. Evaluation and
  fitting refusals propagate; there is no silent low-ka fallback or unloaded
  assembly success. `fs-bem` is not a production dependency of this crate.
  G0 `thin_plate_refuses_instead_of_low_ka_fallback` exercises this refusal.
  The existing compact one-branch load keeps ka<=0.5 and dt*fmax<=0.1;
  its analytic leading resistance and inertance coefficients are tested at
  1e-14, with a 512-frequency grid targeting <2.2% complex and <4.2%
  resistance error versus the closed form (native measurement pending).
  These are sampled discrepancies, not continuum enclosures. Linear plates
  accept in-plane pretension and clamped edges. Von Karman is the
  isotropic SS sine path when `e1 = e2` and edges are free to
  rotate; clamped or orthotropic bending uses DKT-sampled
  displacement with the same sine Airy basis and the complete oriented
  membrane stiffness `A = 12 D / h²`
  (`von_karman_sampled_plate_with_membrane`). Independent shear,
  Poisson and rotated normal/shear coupling survive the reduction.
  The compact far-field observer is
  the baffled on-axis piston `p = ρ A ÿ / (2 π r)`, the same
  half-space as the self-load. A full BEM observer is the
  solver-neutral `DirectionalFarFieldTable` /
  `far_field_observer_pressure` (`p = F e^{ikr}/r`); `fs-bem`
  produces `FarFieldTable` as a **dev-dependency** only (the
  production cycle `couple → bem → solver → feec → couple` is
  refused).
  `CompactBody::drive` reuses the shared exact-ZOH modal transition under
  held force, converting through mass-normalized coordinates. Acceleration,
  volume velocity and pressure describe the completed endpoint. G1/G3 tests
  compare analytical damped/undamped motion and held-force subdivision,
  including a step beyond the former Euler stability limit. A retained piston
  filter is driven by the exact same-step mean surface velocity; its direct
  reaction is solved algebraically before either history is accepted. The
  resulting reaction work is `-dt * face_area * pressure * mean_surface_velocity`,
  so a resistive load cannot inject energy through a lagged velocity. The
  scalar transition and filter feedthrough are prepared once, with mechanical
  cache refresh when the public frequency or damping ratio changes; mass and
  signed modal area are applied at their current values. The mean relation
  retains small increments without subtracting absolute endpoint positions.
  Negative or nonfinite instantaneous filter resistance refuses. This is a
  power-consistent sampled interconnection, not an exact continuous coupled
  solution or a global positive-real certificate for the fitted filter.
  Compact motion and sampled radiation-filter memory are proposed together.
  `drive_and_radiate` admits positive finite observer density/distance and finite
  pressure before publishing either history. The independent `PlateBank` path
  prepares every mode and checks the final pressure sum before publishing all
  modes together; a later mode refusal leaves earlier accepted modes unchanged.
  Refused motion, filtering or observation can therefore be retried from the
  same accepted physical histories. Filter coefficients stay in place, while a
  second mutable filter-state buffer and O(mode-count) scalar proposals are
  reused after initial allocation. This preserves exact scalar motion under
  the solved held force and pressure summation order; it does not promote sampled independent radiation
  loads to the joint continuous-radiation model. Focused tests in
  `thin_plate/bank/tests.rs` cover late-mode refusal, pressure overflow and
  bitwise-equivalent retry, including nonzero retained filter histories.
  Manufactured resistive and first-order passive loads additionally close
  structural work, retained radiation storage and physical heat; a strongly
  resistive negative control demonstrates the former lagged instability.
  The finished pressure history is passed through ISO 9613-1
  absorption (`air_path`) with the assembly's explicit humidity;
  Stokes–Kirchhoff is only the fallback outside the ISO window.
  ISO is applied only with an `fs_material::iso9613::TerrestrialAir`
  admission bound to the gas state. The scenario's legacy moist air makes its
  own. An explicit `GasAbsorptionModel::Iso9613Air` carries one made from the
  gas's composition, so heliox, oxygen-poor or CO2-rich gas refuses with
  `AcousticRealizeError::NonTerrestrialAir`, naming the species.
  The receipt's `relaxation_omitted` is
  `OmissionBound::StandardAirConvention`: standard air's worst ISO absorption.
  It is a convention and not a bound for a gas that relaxes more strongly than
  air (CO2-rich, for example).
  `air_path::omitted_observer_absorption(record, range, band_top)` is the
  figure a render that applies no law (or Stokes-Kirchhoff only) leaves out,
  read from the listener gas's record, never inferred from its state: ISO
  9613-1 for a legacy record or terrestrial-air composition; Stokes-Kirchhoff,
  complete, for monatomic species only; the relaxation law with ISO's rates
  held (declared) for N2 and O2 off ISO's fractions; otherwise `Unbounded`,
  naming the species (heliox's He, CO2 above ISO's ceiling, SF6). The band top
  is `observer_band_top_hz`: the declared band, Nyquist and 20 kHz
  (`OBSERVER_BAND_TOP_CAP_HZ`), whichever is lowest.
  `ObserverOmission::admit` refuses with
  `AcousticRealizeError::ObserverOmission` when the omitted figure exceeds
  `OBSERVER_OMISSION_TOLERANCE_DB` (1.0 dB) or is unbounded, unless the input
  declares a law (`observer_law none` or `stokes-kirchhoff`), which accepts the
  omission with its receipt. In air at 20 °C, 50 % RH and 1 atm, at 20 kHz, a
  1 m listener omits 0.524 dB and passes; a 2 m listener omits 1.048 dB and
  refuses. A gas without a record naming its composition (a resolved field
  cell) is `Unrecorded` and is treated like `Unbounded`.
  `GasAbsorptionModel::MolecularRelaxation` applies the gas's own relaxation
  law (`fs_material::relaxation`). It is evaluated from the gas's composition
  and its species' data, and supplied as `RealizeOptions::observer_relaxation`
  at exactly the realized gas state. It is applied as the causal transfer
  `air_path::propagate_pressure_history_relaxing`: absorption and dispersion
  come from one wavenumber and Stokes-Kirchhoff is counted once. Batch and
  streaming propagation share one energy-bounded FIR; batch uses the exact
  reciprocal-dt clock and linear overlap-add convolution. Finite records are
  zero-extended, so late arrivals cannot wrap into the start, and the tail
  beyond the supplied record is discarded. The returned receiver delay is
  `r/c_0 - FIR_DC_group_delay`, including the FIR's realization latency, and
  remains signed; it is not applied to the source-time history here. Invalid
  input, response or output refuses without changing that history. The receipt
  names the law, with nothing omitted.
  These refuse rather than being ignored:
  - the law without its state;
  - a state of another gas;
  - a state supplied with any other law.
  `xtask/tests/observer_relaxation_amount.rs` gates the amount on a rendered
  plucked string at 100 m. The humidity change is gated against ISO 9613-1,
  and the oxygen change against the law with declared rates.
  The shared `pcm_wav::baffled::BaffledPressure::new_multirate_gas_differentiating`
  composes these same `ObserverFilter` laws on every original surface-to-receiver
  path before signed modal summation. Only bit-equal computed distances share a
  prepared filter; there is no common-distance filter on the summed pressure.
  Each path retains its `rho*dA/(2*pi*r)` gain, source filter receipt and residual
  delay `r/c - filter.low_frequency_delay_s`. A negative residual receives the
  smallest common whole-sample causal padding across ALL paths; that added
  numerical latency is reported separately and preserves relative path timing.
  The existing 65-tap sinc contributes its original 32-sample interpolation
  latency. The source response's own delay is counted once, including dispersion
  and its numerical padding; FIR support is not physical travel time.
  Explicit distinct-path, retained-response and convolution/projection budgets
  bound cold preparation, with cancellation between bounded designs and surface
  tiles. The actual pressure loop and clone retain the original receiver owner.
  The reference cold budgets are 16,384 distinct paths, 4,000,000 source-filter
  taps and 500,000,000 convolution/projection terms. These are admission limits,
  not a measured real-time guarantee. The receiver exposes its actual multiply
  count per sample; the live-piano clock reports that count and the physical and
  numerical delays separately. Larger responses, including some reduced-pressure
  or higher-rate configurations, can exceed the declared limits and refuse;
  the caller may select larger explicit cold budgets in the native API.
  The law is explicit: classical transport does not infer molecular relaxation,
  ISO requires its state-bound terrestrial-air admission, and a supplied
  molecular state must belong to the same gas. Prepared live-piano atmospheres
  now use this path with the generic classical law; native configurations can
  supply another admitted law or retain the lossless comparison receiver.
  This remains fixed-gas, infinite-baffle, one-way observation, with original
  source-law, filter-tail, temporal and spatial approximation limits.

  `air_path::fluid::{FluidObserverFilter, propagate_pressure_history}` uses
  the original material `ClassicalBulkAbsorption` for liquids and real gases,
  reusing the existing response design, FIR stepping and linear convolution.
  The sample clock remains exactly `1/dt`. Numerical tail-energy tolerance is
  explicit and is not derived from terrestrial-air uncertainty. Admission
  of the classical approximation requires `omega*tau <= 0.01` through Nyquist;
  this supported-regime threshold is not an error bound. The retained material,
  selected model, actual FIR DC
  delay, remaining `r/c0 - delay`, numerical padding and tail energy are
  reported. Nonfinite input/output refuses atomically. Existing constructors
  retain the classical finite-band minimum-phase approximation.

  `FluidPropagationModel::ThermoviscousDispersion` is explicitly selected by
  `FluidObserverFilter::with_model` or
  `air_path::fluid::propagate_pressure_history_with_model`. Its complex
  acoustic root supplies exp(-Im(k) r) exp(-i (Re(k)-omega/c0) r) for the FFT
  synthesis convention. The original response owner retains this phase through
  Nyquist alignment, precursor padding and truncation; minimum phase is not
  reconstructed from its magnitude. `r/c0` is a reference delay, not a finite
  wavefront of the parabolic Navier–Stokes–Fourier equations. Branch refusals
  propagate; no molecular spectrum, boundary thermal mode or temperature
  evolution is added. This full acoustic branch admits strong absorption beyond
  the classical `omega*tau <= 0.01` threshold; the maximum parameter remains an
  explicit diagnostic. Every sampled branch must resolve and the original FIR
  energy, length and representability limits still apply. Material frequency
  support and constant-coefficient Newtonian-continuum applicability remain the
  caller's responsibility; this is no rarefied-fluid or experimental-validation
  claim. Independent viscous-root tests distinguish the strong-loss magnitude
  and phase from quadratic attenuation, retain batch/streaming parity and
  atomic acoustic/thermal branch refusal, and keep the classical guard intact.

  `pcm_wav::baffled::BaffledPressure::new_multirate_fluid_differentiating`
  composes that same fluid filter on each original surface path. Gas and
  general-fluid paths share one preparation owner for exact-distance caching,
  original density gain, signed surface projection, common causal padding,
  cancellation, work limits and allocation-free stepping. It retains the
  original infinite-baffle and prescribed-motion assumptions. G1/G3 tests
  compare signed complex path sums, delay closure and actual pressure changes
  against independently calculated bulk-viscosity attenuation; all original
  gas path semantics remain in the same owner. Explicit `_with_model`
  variants on this surface constructor and the circular-opening constructor
  preserve the selected acoustic dispersion on every original path. Rayleigh
  source normalization and boundary assumptions remain unchanged. Six added
  regressions cover independent viscous complex phase, signed surface paths,
  batch/stream timing, exact-clock admission, path composition and legacy
  constructors. These integrated regressions remain unrun because the execution
  environment became unavailable; they are not a new full-NSF boundary solve.

  `circular_openings_multirate_fluid_differentiating` exposes the same liquid
  and real-gas propagation to native volume-flow apertures. Circular quadrature,
  overlap/band/clock admission and signed input ordering share their original
  gas-path owner. `align_fluid_latency` aligns separately prepared receivers on
  the same output/reconstruction clock without changing physical path receipts.
  It refuses stale padding, decreasing latency, history overflow, directional
  velocity banks, and any accepted trace, including a silent multirate packet
  before pressure arrives. Rejected packets remain cold; clones preserve the
  accepted-clock state. The existing cold-only `differentiate_pressure` change
  likewise refuses after stepping. G1/G3 tests cover independent piston
  directivity and signed split-versus-joint transfer and waveform equivalence;
  G0 verifies refused reconfiguration leaves the live history and report intact.

  Offline acoustic assemblies apply every selected observer law through
  `air_path::propagate_pressure_history_with_model`: ISO's complex Debye
  phase, the classical law's declared in-band minimum phase, and molecular
  relaxation all share the streaming response designer and linear batch
  convolution. The batch clock remains the exact reciprocal of dt. Short and
  non-power-of-two histories are filtered, with zero input before the record.
  The signed remaining receiver delay has the same semantics as above; the
  absorption receipt identifies the law and omissions, not an absolute arrival.
  The explicitly magnitude-only `absorb_pressure_history_with_model` API
  remains available for that declared comparison; assembly defaults no longer
  route through it.

  Streaming renderers apply the same laws through `air_path::ObserverFilter`,
  a causal FIR built once per law, range and clock. The reed performance
  (`render::schedule::reed`) applies it to both observers, the legacy
  `listener_m` proxy and the `baffled_outlet` receiver. It chooses as above:
  ISO 9613-1 for a legacy (T, p, RH) medium, Stokes-Kirchhoff (with the
  omission receipted) for a supplied mixture or thermochemical medium, or the
  `observer_law none|iso9613|stokes-kirchhoff|relaxation
  iso9613-terrestrial-air` record. `none` (`GasAbsorptionModel::Lossless`) is
  a receipted lossless choice for comparisons with lossless analytic oracles,
  never a default. `ReedPerformance::observer` returns the receipt.
  The same shared input in native reed and plate performances also accepts
  `observer_law relaxation declared ...`, using
  `fs-material::relaxation::input::read_relaxation`'s original bounded
  species/reservoir/rate grammar. The law evaluates against the exact exterior
  mixture or thermochemical record; missing species/partners, a legacy
  pseudocomponent or a rate outside its temperature domain refuses. Declared
  data remains authored constitutive input, not measured calibration. Each
  open mouth keeps its own range, causal dispersion and receiver delay; no
  second classical or ISO absorption is added. Existing defaults are unchanged.
  How each law is made causal:
  - Molecular relaxation: its own causal transfer.
  - ISO 9613-1: each O2 and N2 term by its exact complex Debye
    (Kramers-Kronig) pair, so the magnitude is ISO's exactly and the dispersion
    is its own. The classical term is the in-band minimum phase.
  - Stokes-Kirchhoff: the in-band minimum phase. The receipt names the band
    where `omega tau_v << 1` fails.
  A minimum phase built from the in-band magnitude is not used where
  dispersion exists, because it cannot see absorption above Nyquist. For air
  at RH 50 % and 30 m, against the direct causal relaxation transfer:

  | Grid Nyquist | Worst in-band complex error | LF delay (physical 34.5 us) |
  |---|---|---|
  | 24 kHz | 2.0 | 20.2 us |
  | 48 kHz | 1.07 | 27.4 us |
  | 192 kHz | 0.11 | 33.8 us |

  A dispersive transfer sampled only to Nyquist rings before `t = 0`. Its
  excess phase at Nyquist is nonzero, so the periodic DFT spectrum jumps there,
  and no zero-latency causal FIR drops the precursor below a fixed fraction.
  That is a property of band-limited dispersive transfers, not of one law. Its
  FIR is therefore delayed by a sub-sample alignment `d` (making the Nyquist
  bin real) plus the fewest whole samples `D` for the precursor to meet the
  certificate. Each receipt's certificate (`budget_energy`) is
  `min(OBSERVER_TRUNCATION_ENERGY = 1e-9, OBSERVER_UNCERTAINTY_MARGIN = 1e-3 x
  alpha_uncertainty_energy)` of the response energy. `alpha_uncertainty_energy`
  is what ISO 9613-1's stated +-10 % on alpha is in the same measure. The
  dropped energy therefore sits three orders below the law's own uncertainty
  at every range, and short ranges pay more padding. The receipt prints `d`,
  `D`, the budget, the latency and the realized dropped energy. There is no
  fixed range floor: a baffled receiver whose own delay cannot absorb the
  filter's refuses (`ObserverLatency`, naming `observer_law none`), wherever the
  physics puts that for its law, gas, clock and range. A 5 cm receiver at
  48 kHz in air refuses.
  Precursor energy left before `-D`, ISO companion, air, 48 kHz (prototype):

  | Range | Alignment | D = 16 | D = 64 | D = 256 | D = 1024 |
  |---|---|---|---|---|---|
  | 0.5 m | none | 3.1e-5 | 8.0e-6 | 1.9e-6 | 3.7e-7 |
  | 0.5 m | d | 1.7e-9 | 3.3e-11 | 5.5e-13 | 9.4e-15 |
  | 30 m | none | 2.7e-4 | 7.1e-5 | 1.7e-5 | 3.4e-6 |
  | 30 m | d | 3.0e-7 | 6.3e-9 | 1.1e-10 | 1.8e-12 |
  The filter is referenced to its front. The baffled receiver's propagation
  delays use the receipt's `delay_speed_m_s`
  (`RegionalDuct::baffled_outlet_delayed_at`), so receiver delay plus the
  filter's low-frequency delay (dispersion and latency) is `r / c_0` exactly
  once. A receiver too close to absorb it refuses with
  `AcousticRealizeError::ObserverLatency` (range, available and needed
  samples). The law is applied at the receiver's distance from the outlet
  centre; `range_spread_m` bounds the difference over the Rayleigh elements.
  The proxy has no propagation delay, so its filter delay is uncompensated.
  `ReedObserverReceipt::consumer_delay_s` is the delay a consumer must apply
  to place the proxy's radiated part as the baffled path is placed: both
  total `r / c_0`. The proxy's sample is bore pressure plus radiation (plates
  and jet at `listener_m`). The law acts on the radiation only
  (`ReedBoreVoice::last_observer_parts`); the bore pressure travels no
  distance and passes unfiltered and undelayed
  (`internal_part_unpropagated`). Under `observer_law none` the proxy is the
  voice's own sum, bit for bit.
  The legacy jet contribution is an explicitly on-axis compact FORCE dipole:
  `F = (p_m-p_bore) width opening` is already in newtons, and the shared
  `fs-aeroac::curle3d::CompactForceObserver` supplies `Fdot/(4 pi c r)` without
  another density multiplier. It is prepared once from the inlet material;
  density still controls the reed flow and bore impedance. The proxy keeps its
  historical far-field/source-clock scope; its consumer supplies propagation.
  Plate volume-acceleration monopoles retain their required density factor.
  The physical baffled-outlet path instead propagates its outgoing duct wave
  through the existing exterior-fluid observer; an enclosed reed's proxy is
  not added as an independently exposed exterior source.
  `ReedPerformance::observer_admission` reports the listener gas's figure
  (`air_path::omitted_observer_absorption`, from its record; a field cell is
  unrecorded) at the observer's range and band top (the receiver's declared
  band, or Nyquist for the proxy, capped at 20 kHz), and what the applied law
  leaves out: nothing under ISO or the relaxation law, the relaxation part
  under Stokes-Kirchhoff, everything under `none`. It refuses an undeclared
  Stokes-Kirchhoff law whose omission is unbounded, unrecorded or above
  1.0 dB. The library builds regardless; `music_render wind` calls it before
  any output exists and writes the receipt as `reed_input.observer`. The five
  reed CO2 fixtures declare `observer_law stokes-kirchhoff`, the law they
  already had by default, so their sound is unchanged.
  Authored string `ζ` stays a constant modal ratio; a relaxation spectrum
  is present only when the caller explicitly selects `PronyBending`.
  An optional `HelmholtzCavity` faces the plate monopoles as a
  flow-driven pHS whose damper is compact-mouth `Re Z_rad(ω₀)`.
- Bow roughness: an optional `ContactTexture` is a declared
  self-affine height spectrum. `fs-tribo::surface_excitation` samples
  it; the height perturbs the normal load that Stribeck sees. No
  measured rosin is invented.
- Three-way clock: when `EA > 0` the string member is Kirchhoff–
  Carrier with contact inside `H`. Polarizations share the plate and
  the duct FIR. The certified three-pHS Dirac interconnection is
  `common_flow_dirac(moving_end_waveguide, transformer(plate, load))`
  in `fs-phs` (free-end attachment port, area transformer, Kirchhoff
  star). The complete plate image retains its Helmholtz cavity, material
  memory and exterior histories before an optional `acoustic_chain` duct
  joins its independent surface port. A cavity and duct can therefore load
  the same body without downgrading the string to one-way bridge forcing.
  Cavity construction uses the same `PlateBank::attach_cavity` owner on all
  paths; a cavity without a structural attachment refuses instead of being
  silently omitted. `PrestressedString.moving_end` selects
  that clock. Bow, obstacle stations, and blow/reed are leftover
  ports of `join_port` (Stribeck and Hunt–Crossley stay port
  forces). Both moving-end ODE and descriptor clocks bind obstacle stations
  through the shared contact admission and retain each physical line measure,
  gap and law. `fs-dcontact` supplies the instantaneous normal reaction and
  non-attractive Hunt–Crossley unloading at the actual conjugate port velocity;
  `fs-tribo` supplies the existing declared scalar friction ramp using the same
  elastic normal capacity as fixed-end contact. Stations are never relocated
  and missing measures or gaps refuse. G1/G3 checks cover force-law parity,
  nonzero loss, station-measure partition invariance and friction's effect on
  observer pressure in both clocks, including retained Prony memory. The held-force clock
  is unchanged; this wiring adds no finite-step contact-work closure claim.
  Von Karman is the same plate pHS with quartic storage.
  `EA > 0` is `kirchhoff_carrier_moving_end` on the same
  free-fixed port. Fixed-fixed sines keep the one-way bridge
  force; with a plate the duct is still the ODE chain
  (bridge force on the leftover plate port), not the FIR line.
  Without a plate the same cylindrical chain is still the ODE
  clock. A cylindrical bore (uniform or stepped) is the
  `acoustic_chain` LC ladder with an all-regime
  `ViscothermalPin` from `GasState` (`μ`, `γ`, `Pr`;
  series `R` and thermal `G` are three-term Foster networks
  collocated to Bessel Zwikker–Kosten `F(r_v)`, all shear
  numbers). A single widening taper is `spherical_cone`
  (`ψ = x p` plus the one inlet Euler near-field shunt).
  Any other taper, alone or in a multi-section chain, is
  frustum LC cells (exact frustum inertance and volume),
  frankensim-d1j6s. A `Closed` termination is the chain's
  rigid outlet (`ChainOutlet::Rigid`, `U = 0`), matching
  the characteristic line's `Termination::Closed`; before
  frankensim-gs17m it was the pressure-release end, an
  octave away.
  Open tone holes are `AcousticTap` side branches on
  the ODE path: a compact neck is one inductor; a long
  neck (`kℓ > 0.2` or `ℓ > 4b`) is the same 2-cell LC
  line plus flanged mouth the TMM chimney already is.
  The same `ViscothermalPin` sits on a lumped neck or
  on each chimney cell of a line.
  A scenario `LocallyReactingWall` is the ODE `WallPin`
  (per-cell shunt LC, `A_w = 2π a dx` with slant on a
  taper) and the TMM `Y' = 2π a slant / Z'` addend on the
  FIR `Z(ω)` / `R(ω)` path, including tone-hole chimneys.
  Same numbers both clocks. The one-pole `TravelingWaveLine`
  is not that path.
  A closed pad is the TMM cavity compliance on the station
  cell; `foster_branches > 0` adds the bore's thermal
  Foster on that remaining `C`. A scenario `ToneHole.open_fraction`
  σ is the same admittance mix
  `Y = σ Y_open + (1−σ) Y_closed` the TMM `HoleState::Vent`
  already is.
  A quasistatic reed is the Bernoulli port on that
  inlet; a massive reed is `mass_spring_damper` plus the same
  jet, face-velocity flow, and a Hunt–Crossley lay as a
  dissipative port force (not a term in H). With a linear plate
  the mouth is a transformer. An unflanged ODE mouth is the
  same compact `(R, X)` as the TMM: `Re Z_rad` plus mass
  `Δℓ = 0.6133 a` on the last flux. A flanged mouth is
  `0.8216 a` (Rayleigh piston above `ka = 1` on the TMM
  path; the unflanged Nyquist ceiling does not apply).
  `foster_branches > 0` replaces the pin `Re Z_rad` with
  a Foster match of that load across a band (unflanged
  samples stay under `ka = 1`; flanged may use the
  piston). Tabulated `R(ω)` remains the FIR path.
  An open chimney with
  `foster_branches > 0` carries the bore's Foster series
  on a lumped neck, or series plus thermal Foster on
  each cell of a long-neck line (extra states). Frequency-by-frequency
  Bessel TMM (`LossModel::Bessel`, spherical-wave cones
  with local-radius lossy substations) is the FIR path. 3-D jet broadband remains a
  no-claim. Plate ×
  Helmholtz in `cavity_phs` is the same transformer, not a
  staggered pair of steps. Its classical and ISO-air observers use the same
  causal complex transfer and linear batch convolution as acoustic assemblies.
  `realize_plate_cavity_with_observer` additionally accepts an explicitly
  evaluated molecular-relaxation state belonging to the supplied gas, or
  lossless observation. These source-time histories include FIR latency;
  remaining receiver delay is not placed on an absolute arrival clock.
  Exterior molecular propagation does not change the finite cavity's existing
  equilibrium lumped gas law.
  `realize_plate_cavity_compact_with_observer` observes both the plate and the
  original neck for a compact enclosure on one infinite rigid baffle. Its
  outward volume flows are summed with the signs of the reciprocal coupling:
  positive plate flow expands the cavity, and the native neck-momentum gradient
  is inward flow. Pressure is `rho d(U_plate + U_neck_out)/dt / (2 pi r)`, the
  compact Rayleigh limit (Mast & Yu 2005, Eq. 1,
  https://doi.org/10.1121/1.2108997). Dimensions and separation must be small
  compared with wavelength and observer distance; unequal path delay and
  directivity are not modeled. An unflanged outer neck refuses this common-baffle
  path. Existing plate-only observers preserve their semantics. Both paths
  prepare the existing exact ZOH owner once and reuse state buffers, including
  its original integrated heat and work checks. Observers do not add plate or
  mutual radiation reaction, change the resonance-pinned neck loss, or promote
  the linear cavity to a finite-amplitude gas model. G1/G3 regressions compare
  independent coupled harmonic equations, low-frequency volume cancellation,
  density/mass similarity, and the previous exact-ZOH plate trajectory.
  Its modal damping inputs are complete: it never selects aluminum/steel from
  density or takes the solid's temperature from the gas. `PlateCavitySpec::bind_plate`
  uses the ordinary `certified_radiators` owner to rebuild frequency, material
  damping, force participation and monopole area together. Force and flow ports
  are normalized by the same square root of modal mass. Source-bound callers
  pass their `ResolvedPlateSpecimen::plate()` and retain that specimen's receipts;
  this numeric adapter does not mint material authority. Thermal and viscous loss,
  including prestress dilution, are included once by the plate owner. Nonlinear
  plates refuse at this linear adapter. Failed admission preserves previous modes.
  G1/G3 tests carry compiled source curves through this pHS pressure path at two
  independent solid temperatures on both sides of the former density threshold,
  compare damping with thermal diffusion, mean-curvature participation and
  Kelvin–Voigt arithmetic, and exercise repeat
  determinism and failed rebinding. They are synthetic implementation checks;
  this adds neither thermal evolution nor exterior plate-radiation back-reaction.
- Von Karman geometric nonlinearity is the isotropic simply-supported
  analytic primitive when those hypotheses hold, and the FE-sampled
  Airy construction otherwise. The full oriented plane-stress membrane
  tensor is retained; the stress basis remains sine Airy (in-plane movable).
- Reed observer adds the compact far-field dipole of the slit force
  `F = Δp · w · y`. That is not 3-D jet broadband and not the
  fs-aeroac 2-D Curle Hankel kernel. 3-D jet broadband remains a
  no-claim. There is no instrument crate.
- `encode_pcm16_wav` maps pascals through a declared full-scale. It never
  peak-normalizes.

## Invariants

- The Dirac interconnection conserves interface power EXACTLY (to roundoff) —
  the G0 law; incompatible ports are refused.
- Vibroacoustic power identities are exact by derivation and asserted
  at solver roundoff; the input-balance identity is TAUTOLOGICAL for
  any consistently assembled system, so the one-sided interface-normal
  mutation is caught by the `interface == cavity` cross-row residual
  (measured alarm 2.0 vs 1e-16). Repeat solves are bitwise identical.
- Every admitted v2 schema has a shape-compatible power contraction, a checked
  watt-dimensional effort/flow product, and energy. Rotational, electrical,
  and chemical PR-2 kinds additionally require angular momentum, electric
  charge, and amount respectively; legacy seed role vectors remain unchanged.
  Field duality records separate effort/flow function-space roles,
  includes its integration-measure dimensions in the checked product, and
  explicitly assigns the measure to the density side before kind-specific
  dimension comparison.
- Stable relation IDs cannot alias their owned port/boundary IDs. A
  conservative pair requires distinct port IDs, matching physical metadata and
  clock/timestamp, plus outward-from-owner conventions on both sides. PR-1
  refuses common-frame orientations until a public pullback can prove them.
- Storage, dissipation, and source ownership remain separate typed claims. A
  conservative junction cannot stand in for any of them; a dissipative
  relation cannot exist without an evidence reference; a source cannot exist
  without a named accounting boundary whose coordinate/sign binding exactly
  matches the port.
- Every admitted stream has a non-empty duplicate-free constituent axis, one
  shared outward sign convention, finite fixed-dimension rates, fixed bundled
  conservation roles, and exactly one energy chart. Alternate charts cannot
  enter without exact scalar identities and context-bound durable evidence;
  the actual rate axis must equal the proof-bound axis, chemical power is owned
  exactly once, and pressure/deviatoric work uses the single normative
  Cauchy-tension-positive integrated outward-power convention.
- Every accounting report is tied to one exact logical interval, subject-bound
  evidence roles, exact initial/final endpoints, one canonical expected-row
  set, one element axis, one charge policy, and one entropy convention.
  Caller order cannot change arithmetic: rows are sorted by stable contribution
  ID and reduced with a fixed-order compensated sum.
- Global first-law, mass, element, and charge closure use two-sided typed
  tolerances. Entropy generation is never absolute-valued or clamped: only
  `S_gen >= -entropy_tolerance` passes. Individual boundary contributions may
  be negative. Internal pair closure is a separate equality gate on every
  audited axis, and violations retain the exact pair, ports, and rows rather
  than guessing a single cause from a global residual.
- Inventory snapshots and filled contributions retain the exact energy
  reference, element schema, charge policy/basis, and entropy convention under
  which they were admitted; audit refuses same-interval context rebinding.
- A constant-temperature thermal power row must carry only heat entropy and
  satisfy `energy = T * integrated_heat_entropy` bit-for-bit. Radiation on a
  thermal power row requires profile-bound external evidence. Profile-bound
  thermal and stream heat/radiation integrals remain external-evidence claims.
- The energy audit reports an interface-balance failure exactly when some
  caller-supplied exchange has absolute imbalance above `tol` or is non-finite.
- On the added-mass fixture (`μ ≥ 1`): naive staggering (`ω = 1`) diverges while
  Aitken-relaxed coupling converges to `x* = c/(1+μ)`; Aitken never takes more
  steps than a stable fixed under-relaxation.

## Error model

DEF-05 modal apertures now admit relative tolerances strictly below 1/128.
Work and flow closure use the actual declared relative ceiling, with no 128x
post-solve allowance. `ModalApertureFrame` records the work and flow
`SolveCertificate`s (joules and m^3/s), including the observed residual, scale,
requested/effective budgets, gamma-n arithmetic floor and the actual scalar
residual-evaluation count. These explicit budgets never auto-raise below their
floor; a refusal retains `fs_phs::PhsError` through
`AcousticRealizeError::AlgebraicSolve`. No state or sample clock is published.
Work's floor uses `gamma_(4*modes+12)` times the absolute energy/work operands;
flow uses `gamma_(2*modes+8)` times its absolute flow operands. These are
conditional arithmetic diagnostics, not bounds on material, flow-model or
time-discretization error. G0/G4 `modal::refusal_tests` tests the 1/128 boundary
and a 1% unaccounted external-force work defect through the production gate.

Structured `CoupleError`; no panics. In addition to legacy incompatible kinds,
errors cover invalid IDs/empty shapes, pointwise or field-measure dimension
overflow, non-power products, kind-specific or shape/pairing mismatch, missing
kind-required conservation roles, schema-only kinds sent through the raw
interconnection oracle, localized schema conjugacy mismatch, identity aliasing,
accounting-boundary coordinate mismatch, scalar/non-scalar misuse, wrong
added-mass port kind, and non-finite schema-bound values. Stream errors localize
non-finite fields/components, empty or duplicate constituent axes, non-outward
orientation, proof-bound/actual axis mismatch, species-potential accounting on
an element axis, chart/proof binding or identity mismatches, non-positive
density, chemical double counting, and declared/chart energy-rate disagreement.
Window errors additionally localize clock/interval/endpoint mismatches,
evidence-role/subject mismatches, malformed or one-sided internal exchanges,
duplicate rows/ports/elements, missing or unexpected manifest rows, foreign
stream context, missing stream/chemical element schema, unproved charge policy,
negative inventory/tolerance, nonpositive temperatures, incompatible
power-port transfer axes, absent temperature evidence, nonfinite integrated
values, and element/charge axis violations. Physical conservation or
second-law failures are report violations, not construction errors.

## Determinism class

Fully deterministic: schema admission canonicalizes roles by enum order,
stream admission canonicalizes constituents by typed ID, and window admission
canonicalizes manifest, per-exchange, port, and element ledgers by stable typed
identity. Exact crosswalks use a pinned operation order and bit equality;
window reductions use a fixed-order compensated sum. Interconnection, relation
construction, audit, and iterations are pure functions of their inputs.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

- `mm_line` module (bead 3ez8g.4.2): the multimodal characteristic-line
  runtime for brass. `MmLineBank` realizes one exact-FIR plane-mode
  reflectance line PER VALVE COMBINATION from the fs-duct multimodal
  authority. Its length comes from the shared `LineLadder`
  convergence ladder (DEF-02, IRP .82), not from a fixed round-trip
  multiple; an executed lesson is that a four-trip tail left the
  lowest peak 92 cents flat from wrap-around aliasing. The ladder
  starts at four round trips and doubles until every combo's
  fourth-quarter tail energy is below `tail_tol` and its first three
  peaks move less than `peak_tol_cents`. It refuses with
  `LadderUnconverged` at the ceiling. The line is terminated into the
  selected mouth load. `MmLineBank::new` and `MatrixFirLine::{new,with_ladder}`
  accept the shared `fs_duct::MouthLoad`. Its `Tabulated` arm requires a
  matching mouth plane/radius/gas and complete positive-frequency support
  for the requested ladder, including its ceiling's lowest bin and Nyquist.
  One bank uses one common table, so every combo must match its declared
  inlet-relative axial position; differing path-length combos refuse.
  Per-combo tables remain outside this interface.
  Every retained positive bin comes from the table; the bank receipt reports
  `table_bins == n_fft/2` and zero fallback bins. DC uses the table's declared
  finite static impedance through the same multimodal bore at the existing
  1e-3 rad/s quasi-static probe, with no table extrapolation, default -1, or
  clipping. Strict-load errors propagate, and plane reflectance is checked
  before FIR realization with the existing 1e-12 rounding allowance.
  Both staircase FIR consumers now sample the authority's native modal
  reflection directly, including analytic, tabulated, full flanged and
  declared-DC loads. They do not pass through a full inlet impedance and
  back: pressure antinodes retain finite reflected waves. The plane bank
  imposes zero higher-mode inlet flow; the matrix line uses its constant
  positive real reference for every incident port. Those are distinct source
  closures. Existing source coverage, passivity and FIR convergence gates
  still apply after this conversion.
  `BrassMouthSequence` and the `LipPerformance` WAV renderer likewise use
  the native modal mouth-transfer map against the same actual outgoing
  pressure history. The output FIR keeps finite mouth motion at inlet flow nodes
  without multiplying a divergent mouth/inlet-flow ratio by vanishing inlet
  flow. One immutable prepared geometry serves each path's entire output
  sweep; analytic, strict table, declared-DC and full flanged loading remain
  live. Output mode/tail/prefix gates, valve-history handling and original
  source mechanics retain their existing meaning. The legacy lip microphone
  receives mode-zero flow with its original uniform-disk receiver. An explicit
  `modal_baffle <relative_kernel_tol> <maximum_kernel_terms>` lip record selects
  the full flanged matrix for both feedback and every mouth-flow FIR, together
  with finite-distance Rayleigh observation of all original signed `psi_n/S`
  traces. It requires a linear bore and `termination flanged`; the same gas
  and actual staircase aperture own the load and receiver. No foreign exterior
  or scalar load is silently replaced. Higher zero-net-flow modes still radiate
  with their spatial phase and direction. The existing bandlimited derivative
  and causal path kernels retain their reported latency and clock limits.
  `baffled::modal` uses Gauss-Legendre radial nodes and periodic angular points;
  declared counts are multiples of four in 4..=128 and 32..=512. Independent
  quarter/half/full axis checks compare every actual modal kernel in L1 against
  the unit-flow scale `rho/(2 pi d_min)`, never against a possibly vanishing
  pressure. Fine basis Gram and band-resolution checks refuse underresolution
  without changing the signed shapes. The explicit total work cap covers all
  five preparations. These are measured quadrature changes, not a continuum
  error enclosure; temporal reconstruction, omitted modes and finite-body
  scattering remain separate limits. Source FIR tolerances are never relaxed
  automatically. The example `modal-baffled-lip.performance` is authored
  development geometry and mechanics, not an instrument-realism claim.
  `MouthLoad::Matched` refuses explicitly on these staircase time-domain lines.
  The existing FIR Nyquist rolloff and precursor folding remain in effect;
  full-band source coverage is not a new high-band realization-error claim.
  `MmLoad::TabulatedWithFallback` remains an explicit diagnostic policy.
  A spliced analytic fallback outside the
  table's support is a spectral jump, so the ladder refuses it by
  name. Passivity (`|R| <= 1 + 1e-12` on the 4x-oversampled grid) is
  checked on the stored taps without clipping and refuses otherwise.
  A valve switch rebuilds from those taps. Valve
  switches at block boundaries use the carry-outgoing-history lift
  (in-flight waves persist; the new taps govern future reflections)
  plus an optional linear CROSSFADE from the old line — an
  instantaneous reflectance swap is genuinely discontinuous physics,
  and the crossfade is the no-click lane; every switch is logged with
  its lift name and carried sample count. `MatrixFirLine` is the
  full N^2-FIR bake-off arm with per-mode incoming-energy logging.
  It is sized by the same ladder:
  - DC is the matrix's own low-frequency limit;
  - the convergence tail is the Frobenius sum over all N^2 entries;
  - peak shift and passivity are checked on the plane entry;
  - its `MatrixLineReceipt` records each entry's tail.
  `cup_junction_boundary` + `CupState` is the D18 mouthpiece-cup lumped
  compliance, TRAPEZOIDAL so the runtime matches the analytic shunt
  `Z/(1 + i omega C Z)` to `O((omega dt)^2)`. Realization doctrine
  (the recorded bake-off): at a brass input plane the higher modes
  are deeply evanescent, so the MM-derived `R_00` — which already
  contains every interior mode-conversion round trip — is what a
  plane-driven source sees; the receipt (ml-002) measured the matrix
  arm within 1 cent of the dominant line with a plane-mode energy
  fraction above 0.999 at 9x the cost.

- `brass_loop` module (bead 3ez8g.4.3): the composed brass playing
  loop — no `Trumpet` type. `BrassVoice` = a 1-DOF
  `fs_phs::mass_spring_damper` lip island (outward-striking; parameters
  at reduce-lab card scale with provenance labels on every input) x the
  per-valve multimodal characteristic lines x the D18 cup junction.
  PITCH IS NEVER ASSIGNED: the control surface is blowing pressure,
  lip tension, valve combination — structurally enforced by the
  battery's source grep. THE NUMERICS LESSON (executed): an explicit
  Bernoulli flow computed from the previous sample's bore pressure
  carries a one-sample delay whose phase flips the aperture's
  small-signal conductance ANTI-dissipative above fs/4 — a crook combo
  locked onto a 7.3 kHz parasitic scream; the fix is the closed-form
  IMPLICIT flow/junction solve (the linear cup relation substituted
  into the Bernoulli law gives a quadratic in the flow), leaving no
  delay in the flow path while the lip gap stays explicit (lip mass
  filters it). Passivity on the line taps is additionally enforced on
  a 4x-oversampled DTFT grid (inter-bin Gibbs overshoot makes a
  grid-only |R| <= 1 claim insufficient). Per-block diagnostics carry
  the lock estimate, gap statistics, blow/bore work, and the lip
  island's supply defect.

- `brass_loop::nonlinear` (feature `nonlinear-brass`, bead frankensim-giuhr,
  IRP M2.3-W2a): the same lip island on the shared nonlinear Euler bore
  (`fs_duct::nonlinear::bore`). By default the mouth is a held stagnation
  reservoir at room temperature and room pressure plus the blow. Forward flow leaves it through
  the source's isentropic nozzle at the face pressure; reverse flow leaves the
  bore's static end state into the mouth. The face pressure solves continuity
  between the metered lip flow and the bore's own end wave (a strictly
  monotone residual, bracketed between the closed-end pressure and the mouth
  pressure, Illinois false position to 1e-12 of the bracket). If an auxiliary
  high-pressure trial exhausts the inflow's cold source path, bounded bisection
  first finds a supported upper residual; it never clamps the temperature or
  refuses a supported root merely because that auxiliary trial was too cold.
  Inflow energy is mass flux times the mouth's stagnation enthalpy. The
  mixed-out face's static temperature solves
  `h(T_mouth)-h(T_face)=u_face^2/2` through the shared gas owner's original
  variable-Cp caloric inversion. Continuity and final mass, momentum and energy
  fluxes use that same source-admitted face state; the face's density and sound
  speed reflect its kinetic cooling. Regressions check an independent
  variable-Cp enthalpy integral at elevated pressure, nozzle continuity,
  pressure-wave compatibility, the linear limit, reverse flow and supported
  choked inflow behind an unsupported auxiliary pressure trial. The lip
  advances once per gas step by its exact zero-order-hold split under the face
  pressure. The bell flow is averaged per mechanics sample and observed by the
  multirate baffled receiver.
  - **Finite mouth:** `with_finite_mouth(FiniteMouthSpec)` attaches an authored
    rigid, adiabatic, acoustically compact upstream `fs_gas::Plenum`. Its
    initial absolute pressure and temperature use the same frozen source as
    the bore. The blowing-pressure control sets a held feed reservoir at room
    temperature through an authored effective nozzle area; it never resets
    the stored gas. Zero feed area makes a sealed precharged supply. Actual
    stored pressure drives the lip. Forward and reverse feed use the donor's
    source-owned nozzle and stagnation enthalpy; entering kinetic energy
    thermalizes without being counted twice. Signed lip mass and total energy
    use the identical boundary fluxes and SSP-RK2 stages as the actual bore.
    Repeated stage-zero trials restart from accepted storage. Both candidates
    pass original caloric admission and a compact-volume face-speed CFL of
    0.45 before acceptance. The finite inventory and lip island publish only
    after a successful bore step. The bore's existing post-gas wall-memory
    refusal can still leave a partial step and requires ending the run.
    `FiniteMouthState` exposes pressure, temperature, mass and sensible internal
    energy; `FiniteMouthReport` retains accepted signed feed and native bore
    exchanges, both storage residuals and the maximum stage CFL. Physical
    energy adds the original source reference times mass. Candidate storage
    needs no per-step heap allocation. Regressions exercise an independent
    variable-Cp adiabatic tank law under step refinement, conserved bore-plus-
    storage exchanges, reverse-flow filling/heating, feed reversal, rejected
    and repeated stages, and actual lip forcing without control resets.
    This is not an anatomical lung model: no feed-neck inertia, wall heat or
    moving storage volume is implied. Existing lip swept-volume work remains
    outside the approximation, so these gas budgets do not assert whole-loop
    passivity. Jet/constriction series noise with finite storage refuses in
    either attachment order until its work is closed; source-domain exits
    refuse rather than clamp. Attach once before any gas step.
  - **Lip file:** `bore nonlinear <max_cell_m> <mechanics_per_output<=16>
    <gas_steps_per_mechanics<=64> <wall_poles> <band_lo_hz> <band_hi_hz>
    <wall_error<1>`. It needs `ambient thermochemical ...` with its source packs
    (the bore evolves the ambient gas's own NASA-9 laws) and `termination
    flanged` (its bell is a baffled piston). The feature is off by default:
    a default `music_render` refuses the record by name (tested in a
    default-feature build) and never falls back to the linear bore. Build:
    `cargo build --release -p fs-couple --bin music_render --features
    nonlinear-brass`.
    The optional row `mouth_reservoir <volume_m3> <feed_area_m2>
    <initial_pressure_pa> <initial_temperature_k>` realizes the finite mouth
    above in this same engine. Its initial pressure is absolute; the existing
    blowing-pressure track controls the upstream held feed relative to room.
    Volume, pressure and temperature must be finite and positive; feed area
    is finite and nonnegative. The row requires `bore nonlinear` and refuses
    duplicate records or either series-noise row. The initial temperature is
    admitted by the original source at realization. Without the row the held
    mouth remains the default. `LipDescription::mouth_reservoir` retains the
    authored record; `LipPerformance::finite_mouth_state` and
    `finite_mouth_report` expose the actual accepted owner when the feature
    is built. Parser/render regressions check scalar and feature refusal,
    source admission, sealed-storage independence from feed gestures, and
    positive-area feed filling in the real nonlinear renderer.
    The wall memory is admitted per radius by `fs_phs::ZkDiffusiveBasis` on the
    declared band: on 20-8000 Hz at 1 %, 24 poles resolve radii up to about
    6 mm and 64 (the fs-phs cap) up to about 35 mm; a 61 mm trumpet bell needs
    a declared error of 7e-2 at 64 poles.
  - **Limits:** reverse flow uses the static end state as donor; the cup is a
    lumped volume on the first cell; valves refuse. The field keeps to the
    NASA-9 region that holds the
    room for every species (the McBride air packs split at 1000 K and state no
    caloric continuity across it); a state that would leave it refuses. A test
    builds that air from the data/matdb/seed-v1 McBride tables: the old fixed
    200-1000 K domain refuses it, the room's region admits it, and a room state
    round-trips at 343.24 m/s.
  - **Linear limit** (GO condition, mail 9115): a 1e-8 m^3/s pulse into the
    fixture bore reproduces the multimodal input impedance the linear
    `MmLineBank` is realized from (Bessel loss, flanged end). The first six
    peaks: within 10 cents, 1.5 dB and 10 % in Q (tolerances fixed before the
    run). Measured with 64 poles on 20-8000 Hz (wall
    error 4.9e-3 complex, 6.0e-3 loss): -0.82 to +2.12 cents, -0.10 to -0.25 dB,
    Q -0.6 to -1.8 %.
  - **Cost** (ts1 at load ~100, trumpet Bb4, 4 s at 32 kHz): the linear loop
    takes 132-135 s per audio second, the nonlinear bore about 7000 s, about
    80 % of it in the four-species NASA-9 mixture caloric law (the wall update
    is 8 %). The next lever, fewer caloric-law evaluations per gas step, is the
    fs-gas/fs-thermochem owner's call and is not done. The fixed-lip-input
    trumpet A/B against UIowa is recorded on bead frankensim-giuhr; the graded
    trumpet cards stay on the linear bore until that A/B shows the nonlinear
    bore helps on the origin grader.

- `flue_loop` module (bead 3ez8g.10.3): the composed flue playing loop
  — no `Flute` type. `FlueVoice` = the jet-drive island x the wind
  line. NO CARD -> REFUSE: the island is parameterized only from a jet
  card minted by fs-aeroac's lab (`JetIsland::from_card`), and only a
  card whose claim class can drive a tone is admitted (the tonal
  interim card; the broadband-refusal-boundary card refuses with
  `CardClassCannotDrive`). Card-backed numbers: convection ratio
  `c_p/U = St (h/delta) / (stage + 1/4)` from the stage-I lock, the
  profile ratio `theta/b`, the Reynolds validity band, the claim class;
  authored and labeled: the spatial deflection gain, the onset seed,
  the legacy mouth radiation resistance (unless physical reaction is selected).
  PITCH IS NEVER ASSIGNED: blowing
  pressure and geometry gestures (cut-up, labium offset, jet angle)
  are the controls. Model: Bernoulli jet speed; receptivity = the
  leaky INTEGRAL of the delayed mouth velocity (its 1/omega cancels the
  drive derivative's omega, so the pipe's resonance alone selects the
  mode); the sinuous mode's amplification BAND from the card's theta
  (two low-pass sections at Michalke's tanh neutral `omega_c = U/(4
  theta)` and a high-pass at `omega_c/2` because growth vanishes at
  low Strouhal); labium cut `Q_j = Q_max (1 + tanh((eta - y0)/b))/2`;
  series drive source `rho (Delta_m/S_m) dQ_j/dt` in an implicit mouth
  loop with the mouth inertance, radiation resistance, and the line;
  an EXACT per-block energy ledger (source = stored + numerical +
  radiation + pipe; closure ~1e-15). THREE EXECUTED TRAPS on the way
  to a note: (1) a frequency-flat deflection with the omega of the
  drive derivative locks the highest stage the delay allows (9.4 kHz);
  (2) a continuous onset seed at 1e-3 U, integrated, flaps the jet by
  its own half-width and the labium slope averages to nothing (no
  lock; the seed-independence falsifier in the fixture is the guard:
  a limit cycle's RMS does not move when the seed drops tenfold);
  (3) an amplification band with no low edge lets the integral gain
  win at 6 Hz where the open pipe has no impedance (a relaxation cycle,
  not a note). Historical observations on the flute-class fixture (6 mm bore, 12 kHz
  line for the radiation law's `ka < 1`, 10 x 1 mm flue, 4 mm cut-up,
  tonal interim card): 10-50 Pa lock the fundamental (0.90-1.01 f1),
  110-240 Pa the octave (1.85-2.02 f1), 360 Pa the twelfth (2.84 f1);
  a 5 mm cut-up locks one mode below a 3 mm one at equal pressure; four
  rank voices are bitwise their solo selves in any stepping order;
  replay is bitwise and keyed by Philox voice identity. Those register-pressure
  observations predate the retarded-time correction below and do not qualify
  its corrected feedback phase; they require a fresh sweep. The card's
  Reynolds band (144-264) does not cover playing pressures (Re ~1500),
  so `claim_check` refuses every claim from this fixture: that boundary
  is inherited by the flue-gates bead, not papered over. This chain is
  the DEFAULT, and its gain is AUTHORED: the cards' gain 90 is about
  7-18x the forced jet's response at each note's speaking mode (mail
  7491). It stands until the forced jet passes fg_003, fg_004, fj_001 and
  gas_direction_amount (ruling 7518, mail 7859, bead pexfp.135.1).

  OPT-IN MOMENTUM DRIVE (`FlueVoice::with_momentum_drive`, card row
  `momentum_drive KAPPA` right after `jet_response`): the forced jet's cut flow
  also drives the pipe through the jet's momentum flux, Thwaites and Fletcher
  (1983, JASA 74:400) Eq. (1) p. 401: a series pressure `kappa rho v0 (Q_j -
  Q_rest) / S_p` beside the volume drive `source_coeff dQ_j/dt`. `v0` is the
  bulk jet speed of `q_max = v0 h W`, `S_p` the bore area, `Q_rest` the cut
  flow of the quiescent jet (no source from the steady jet); the term is in the
  same-sample solve and in the `feedback < 1` guard (never loosened). KAPPA is
  required and only `1` (Eq. (1), as printed) or `2.6` (their Sec. II
  pulsating-jet fit, one rig, authors unsure of its meaning) is accepted;
  other numerical weights are refused. A card without the row is bit-identical. Transfer
  caveats (turbulent jet, cut-up/slit 20-40, versus this laminar Bickley jet;
  their end correction 19-30 cm referred to the pipe versus ours of a few
  millimetres) are in the module doc and the provenance sidecar.

  PROFILE MOMENTUM (`FlueVoice::with_bickley_momentum_drive`, card row
  `momentum_drive bickley-profile`) integrates the original station's
  `u=U sech²((y-eta)/b)` instead of multiplying its cut volume by one speed.
  Both volume and momentum use that forced jet's actual station width,
  centre speed and span, including its named spreading law; this selection
  also enables the existing `jet_cut_width station` path. With
  `s=tanh((eta-offset)/b)`, the cut is `Q=U b H (1+s)` and
  `J=rho U² b H (2/3+s-s³/3)`. Full-jet momentum is `4 rho U² b H/3`.
  The original series port takes `(J-J_rest)/S_p`, including its signed work;
  the quiescent cut at the same instantaneous speed/offset supplies J_rest.
  A factored cubic difference retains small cut changes. Its derivative
  `dJ/dQ=rho U (1-s²)` preserves the existing same-sample solver's upper
  feedback bound `rho U/S_p`. There is no fitted multiplier, time filter or
  extra lag. The printed lumped and profile sources exclude one another.
  This assumes deposition of selected streamwise momentum into the bore;
  labium recoil, downstream jet dynamics and experimental sound accuracy
  remain outside the claim. G0/G1/G3 kernel tests compare independent
  velocity-squared quadrature, the analytic derivative, cut symmetry, zero
  deviation and momentum-preserving spreading. Native voice/card tests cover
  shared profile selection, same-sample feedback and source-work accounting.

  OPT-IN FORCED JET (`FlueVoice::with_forced_jet`, card row `jet_response
  forced-bickley SPACING WINDOW`; rulings 7376/7518): the forced Bickley
  jet's complete causal response replaces the band chain and the authored
  gain. The labium displacement is `eta_n = sum_k g_k v_(n-k)` of the exit
  transverse velocity (mouth velocity plus onset seed). `g_k` are the hat
  weights of one step response integrated by `fs_aeroac::forced_jet` at
  the station `cut-up / b` and the Reynolds number of the held pressure
  (the gesture's highest blowing pressure, including a declared source's
  pressure path). The exit is the Poiseuille slit (`b = 0.4 h`, `U0 =
  1.5 U`) behind a slit windway attached first. Behind a developing windway
  attached first it is that exit at the held pressure: `U0` the inviscid
  core `sqrt(2P/rho)` and `b = 0.75 (h - 2 delta* (1 + 1/H))`, the momentum
  flux of the core and its two parabolic wall layers (`H = 2.5`), or the
  Poiseuille exit where the layers have met; a held pressure outside that
  law's validity is refused. This is Verge's matching as Segoufin states it
  (thesis 2000, sec. 5.2): Bernoulli centre speed, momentum conserved between
  exit and jet. Without a windway it is a uniform lip jet (`b = 0.75 h`,
  `U0 = U`). At variable speed, accepted velocity
  intervals retain their own advances of `tau = integral U0/b dt`; a new
  speed never retroactively retimes the history. The explicit fixed-Re
  extension is `dz/dt = (U0(t)/b) A z + B v`, with held interval speed and
  linear physical-time velocity. Local derivative-knot integration avoids
  cancellation of accumulated response integrals. The current velocity's
  affine coefficient still enters the same-sample mouth solve; the response
  is held after its window. Constant positive speed from rest retains the
  old cached hat-weight arithmetic. Zero speed disables input/output and
  freezes this approximate clock, with no claim of molecular decay at rest.
  Trial preparation never commits the clock or velocity; accepted publication
  alone appends an interval. The declared minimum positive speed bounds
  the same number of retained samples (three scalars per interval).
  G0/G1/G3 source tests in `flue_loop::forced_jet::tests` cover original
  uniform arithmetic at three clocks, independent response-integral
  correspondence after speed changes/ring wrap, the instantaneous-retiming
  falsifier, and transactional rest/refusal behavior.
  Register physics is graded on the
  production row 0.1 b and 90 b/U0 (`PRODUCTION_SPEC`). Numerical fixtures
  declare `COARSE_SPEC` 0.25/40, which renders the whistle D5 card within
  0.3 cents and 0.3 dB of the production row. A process memo retains 32
  integrated responses keyed by their exact inputs; the integration is
  deterministic, so reuse is bit-identical.
  REFUSALS: current-sample feedback beyond the solve's uniqueness bound
  asks for refinement (a jet that crosses the cut-up and grows within one
  sample; at 48 kHz over flue 0.1-1 mm, cut-up 0.06-0.8 mm and
  100-4000 Pa, only the 5.3 b and 10.7 b stations at 4 kPa reach it); a
  cut-up shorter than one grid cell of the response refuses at attach; a
  retained history at the declared minimum speed beyond the declared budget
  refuses; so do a non-positive reference pressure or one whose jet is below
  that minimum, the unsteady windway
  (its jet starts from rest, where this frozen-time response needs
  unbounded history), a windway attached after it, and the Nolle/Price
  transfer beside it.
  NOT CLAIMED: the base jet's own spreading over the cut-up, the labium's
  hydrodynamic back-reaction, the Reynolds number of attack and release
  (frozen at the held pressure), and memory beyond the window.
  Of the four gates of the switch, its fg_004 twin
  (`fg_004_forced_jet_*`) passes unchanged, but it LANDS HIGH on the three
  register gates. Their ignored forced-jet twins
  (`fj_001_forced_jet_*`, `fg_003_forced_jet_*`,
  `gas_direction_amount_forced_jet.rs`) carry the measured registers. The
  default switches only when all four pass with their assertions unchanged.

  OPT-IN JET SPREADING (card row `jet_spreading bickley|fletcher-1976`,
  after `jet_response` and after `momentum_drive` when both are present;
  `FlueVoice::with_forced_jet_spreading`, `forced_jet::JetSpreadingLaw`;
  MB 97/125/216, bead pexfp.135.1): the forced jet's response is integrated
  over a widening base (`fs_aeroac::forced_jet::BaseSpreading`), held past
  the station, in a band of `BAND_HALF_HEIGHT` station half-widths rounded
  up to whole spacings. The voice's tanh cut takes the station's half-width,
  `beta(X)` times its own. Its swing `q_max` takes the station's volume
  flux, `beta(X)^(1/2)` times its own. The momentum drive's `v0` is the
  station's centre speed `U / beta(X)^(1/2)`, and the `feedback < 1` guard
  reads that same speed. `bickley`: the laminar Bickley jet at the
  response's Reynolds number, matched at the exit, `beta = (1 + 6 x /
  Re)^(2/3)`. `fletcher-1976`: an organ-flue Pitot measurement (E4 pipe,
  flue 0.25-1.0 mm) carried to the card's jet, not measured on it, Fletcher 1976 eqs. 4-6,
  `phi` 5 deg (Table I), `b = W`, `db/dx = 2 tan phi`. The voice
  provenance and the render sidecar (`flue_input.jet_response.spreading`:
  law, `source`, station width and centre speed, virtual origin, band, cut
  factors) print the law and `JetSpreadingLaw::source()`. On the D4 flute
  card at the coarse row: Bickley `b(X)/b0` 1.0267 and `u_c` 0.9869 `U0`;
  Fletcher 2.5485 and 0.6264 `U0`. Graded through MB (216): Bickley is the
  flute's opt-in. Fletcher sharpens D4 by 13.6 cents and overblows G4 to
  the octave, so it was tried and is not posted. Without the row every
  result is bit-identical
  (`g0_forced_jet_without_a_spreading_row_is_bit_identical`; at 9c04f0065,
  before the jet_noise band of 4ae3c06db, the posted flute D4 and pan-pipe
  D3 mk1 cards rendered byte-identical on this code). REFUSALS: the
  row without `jet_response` before it, `momentum_drive` after it, a law
  other than the two named, trailing fields, and a station half-width that
  is not finite. With the row, the forced jet's first NOT CLAIMED item (the
  base jet's own spreading over the cut-up) is this quasi-parallel law;
  without it, that item stands. NOT CLAIMED: the base's transverse velocity
  and streamwise derivatives (quasi-parallel), the jet past the labium, and
  any spreading law measured on a flute lip jet.

  OPT-IN JET FORCING AREA (card row `jet_forcing_area open-strip`, right
  after a `window embouchure` row, with `jet_response` on the card;
  `FlueVoice::with_jet_forcing_area`; MB 1089/1091/1407, bead pexfp.135.1):
  `v = q / S_m` takes flue width by cut-up as the opening the mouth flow
  crosses. On a recorder that is the window. A flute's flow crosses the
  embouchure's open strip, ALONG by OPEN, which the flue row does not know,
  so the row takes the jet's transverse velocity over that strip. The area
  comes from the window row and is not fitted. Only the jet's velocity moves:
  the source coefficient, the mouth contraction and any authored linear mouth
  loss stay on `S_m`, though the same flow crosses the same strip, so this is
  a partial change, and the render sidecar says so
  (`flue_input.jet_forcing_area`: law, area, mouth area,
  `still_on_mouth_area`). The series mass needs no move: the window row this
  row requires already replaced it. On the held-room flute cards the
  area is 6.44e-5 m^2 against `S_m` 3.68e-5 on G4 and 6.92e-5 against 4.64e-5
  on D4. There the labium swing's fundamental falls from 3.47 to 2.02 h on G4
  and from 3.25 to 2.24 h on D4, and the time past 2 `b_cut` falls from 0.878
  to 0.568 and from 0.898 to 0.836, with both notes in the same register. The
  swing is now below the 3.3-3.95 h de la Cuadra (2005 sec. 2.3) measured on a
  free jet with no labium at Re 200-500. The cards run near Re 1750, so that
  band neither supports nor rules out the row. Graded through MB (1407): G4 3
  to 4 of 14, D4 6 to 5 of 14, because D4's HNR leaves the ruler (34.9 against
  28.4 dB). The spectrum distance falls from 21.05 to 16.24 on G4 and from
  20.04 to 16.25 on D4. Without the row every result is bit-identical
  (`g0_jet_forcing_area_at_the_mouth_area_is_bit_identical`; on this code at
  674860b2e the held-room G4 card without the row rendered wav 0820a975 and
  flows 866aa73f, the posted control's bytes). REFUSALS: the row with no
  `window embouchure` before it, a law other than `open-strip`, trailing
  fields, a card without the forced jet, a late attachment, an area that is
  not finite and positive, and, while rendering, a mouth flow whose speed over
  the named area reaches the gas's speed of sound (MB 1833: PF-661 slide G4 and
  C5 cards with the source term and `mouth_loss none` passed it and returned
  rc 0). The check only reads and refuses, so a render that stays below it is
  unchanged. NOT CLAIMED: a measured flute mouth flow, a uniform velocity over
  the strip, and the terms left on `S_m`.

  OPT-IN JET SOURCE ON THE STRIP (term `source` on the same row,
  `jet_forcing_area open-strip source`; `FlueVoice::with_jet_source_area`; MB
  1497/1783, bead pexfp.135.1): the volume source is driven by the flow the
  jet puts into the mouth, and on a flute that flow crosses the open strip, so
  the term evaluates the default coefficient `rho 0.61 sqrt(A/pi) / A` on the
  strip instead of `S_m`. The coefficient falls by `sqrt(S_m / A)`: x 0.819,
  0.756 and 0.707 on the PF-661 slide D4, G4 and C5 cards. Only the area can be
  argued from geometry. `0.61 sqrt(A/pi)` is the same unflanged end-correction
  law on another opening, not a length measured for a strip. The mouth
  contraction and any authored linear mouth loss stay on `S_m`. The sidecar
  says so: `flue_input.jet_forcing_area` drops the source from
  `still_on_mouth_area` and adds `source_coefficient_on_strip_pa_s2_m3`, and
  `window.jet_source_coefficient` reads `moved to the open strip`. Graded
  through MB (1783): 5, 5 and 9 of 14 on D4, G4 and C5 against 5, 5 and 8 on
  the slide cards, with spectrum, harmonics, envelope and pitch closer on all
  three and centroid inside on D4 and C5; flux leaves the ruler on D4. Moving
  the contraction to the strip as well was rendered on scratch code and graded
  5, 5 and 4, rejected, so it is no term here. A row without the term is
  bit-identical (`g0_jet_source_area_at_the_mouth_area_is_bit_identical` for
  the law at `S_m`). On this code at d2fb1d874 on ts2 the slide G4 card, which
  has the bare row, rendered wav 62caaebe and flows 234a4215, that card's
  bytes, and the D4 card with the term rendered wav bee42500 and flows
  8fee179b, the graded rung's bytes, which were rendered on ts1.
  REFUSALS: a term other than `source`, trailing fields after it, and the
  row's own refusals above. NOT CLAIMED: a measured source length for an
  embouchure, a uniform flow over the strip, and the terms left on `S_m`.

  OPT-IN JET SWING LIMIT (card row `jet_swing_limit fletcher-1976`, after
  `jet_spreading`; `FlueVoice::with_jet_swing_limit`, `jet_swing_limit`; MB
  1784/1816): the forced jet's displacement is linear, and on the organ B4
  card it swings 25 mm across a 7 mm cut-up. The row limits it by Fletcher's
  envelope, his model choice (JASA 60, 1976, pp. 927, 933), carried to this
  jet and not measured on it: exponential growth up to the spread width `W`
  at the lip, then linear growth with an envelope of semiangle about 15 deg;
  he says varying these "has very little effect". `E = W + (X - x_t) tan 15
  deg`, `x_t = X acosh(W w/v) / acosh(y w/v)` from his eq. 9-12, 0 when `W w/v
  <= 1`. `W` is the station's `b(X)`, `X` the cut-up; no coefficient. THE TIME
  FORM IS OURS, NOT PRINTED: one gain `g = min(1, E/y)` on the whole linear
  displacement, stepped at its sign changes, with `y` and `v` its peaks over
  the last two half cycles and `w` from their length; `mu` from the kernel's
  gain `y w/v` read through his cosh, though the kernel's response is not a
  cosh. A saturation of the displacement leaves it crossing zero at the same
  speed and doesn't change the cut's harmonics (open loop on the organ B4
  render). The current-feedback guard still reads the linear slope, which
  bounds the limited one. The render sidecar prints the law, what is printed
  and what is ours (`flue_input.jet_response.swing_limit`). Measured anchors,
  neither on this pipe: de la Cuadra and Fabre, ISMA 2004, amplitude at
  transition about `h` at `Vac/Ub` about 1e-2 and under, and the B4 card's
  `u/U_j` 0.287 is about 30 times outside that; Verge 1995 ch. 4, linear to
  `h/2`. Without the row every result is bit-identical
  (`g0_jet_swing_limit_under_the_spread_width_is_bit_identical` shows the row
  itself changes nothing until the gain leaves 1). REFUSALS: no forced jet, no
  named spreading law, a late attachment, a law other than `fletcher-1976`,
  trailing fields, and a width, cut-up or step not finite and positive. NOT
  CLAIMED: a measured limit on any pipe here, the vortex street Fletcher sets
  aside, and any phase change past the transition.

  OPT-IN JET CUT WIDTH (card row `jet_cut_width station`, after
  `jet_swing_limit`; `FlueVoice::with_jet_cut_width_station`; MB 1784/1816):
  the tanh cut's half-width is half the flue height times `beta(X)`, the width
  at which `2 b H U` equals the exit's bulk flux `U h H`. The prints have one
  `b`, the jet's own sech^2 width at the lip, in the profile, the tanh and the
  swing (Fletcher and Douglas 1980 eq. 1-2, Verge 1995 eq. 4.25, Auvray et al.
  2014), and the kernel already carries it as `b0`. The row gives the cut
  `b0 beta(X)` and the swing `2 b H U` at the station's centre speed. The slope
  at the crossing is the same, so only the saturation differs. `momentum_drive`
  keeps the station's centre speed. Without the row every result is
  bit-identical. REFUSALS: no forced jet, a late attachment, a law other than
  `station`, trailing fields, and `jet_swing_limit` after it.

  JET SPEED EDGE (ruling 8132, bead pexfp.135.2): the declared minimum jet
  speed U_min (`jet_limits`, or 1 m/s on the legacy TMM image) sizes the
  retained history, so every lookup at or above it is retained. A slower
  positive jet is admitted and its delay is never clamped. In the onset,
  before the jet first reaches U_min, a lookup beyond the history that
  precedes the voice's first motion reads the rest state. That is exactly
  zero at any delay: unbounded history at zero memory, bit for bit. The
  forced jet's response reads input older than its first accepted interval
  the same way. Its intervals keep their own advances, so a release usually
  still reaches the window; once a jet below U_min fills the ring with
  intervals that do not, the forced clock freezes as at rest and the voice
  takes the declared state. Any other lookup beyond the history below
  U_min (a release, a dip mid-note, an onset slower than its own transit)
  takes the declared state: no deflection, so the windway's steady flow
  splits at zero deflection, `Q = q_max/2 (1 + tanh(-y0/b))`, and its ramp
  source stays.
  `FlueBlockDiag::below_speed_edge_samples` counts those samples. U_min is
  a declared validity edge, not a stability limit: below it the jet carries
  at most U_min/U of the held flow, and the parallel Bickley jet stays
  unstable down to Re ~ 4 (Tatsumi & Kakutani 1958). Zero blowing keeps its
  own path; the Nolle/Price transfer and the unsteady windway keep their
  own histories.

  Quasi-steady convection reads mouth velocity at the actual retarded time
  `t_n - tau`: the newest committed history entry is `t_(n-1)`, not `t_n`.
  For transit shorter than one sample, the current-velocity interpolation
  term is solved simultaneously with the jet cut-flow and acoustic mouth,
  including any selected contraction loss. Trial evaluations preserve the
  filter memories and consume no extra onset-noise draws. A sufficient global
  uniqueness bound is required before committing the mouth wave; unresolved
  feedback refuses rather than acquiring an artificial extra sample of delay.
  The interpolation regression establishes second-order retarded interpolation,
  not second-order accuracy of the complete jet integrator or measured flute
  fidelity. Gas and windway material changes cross the one-sample boundary
  through the same physical transit calculation.

  Optional `with_viscous_windway(length_m)` binds the same material-gas
  `ViscousSlit` used by reeds to the actual flue exit width/height. Its solved
  bulk speed reaches volume injection, Reynolds number, convective delay,
  amplification and the physical outlet receiver. Schedule admission and
  accepted pressure controls check that actual speed against history/Mach
  bounds; cold attachment also rechecks a previously held pressure. Omission
  preserves the ideal jet. `flue_windway` tests cover viscosity-only material
  swaps at unchanged density/sound speed, actual microphone pressure, control
  refusals and callback partition invariance with simultaneous mouth loss.
  The supplied rectangular channel is a quasi-steady parallel-plate reduction,
  retaining the card profile; exit kinetic head is not acoustic mouth loss.
  No entrance/profile development, windway inertia, turbulent jet, transverse
  flute lip geometry or measured instrument fidelity follows from this option.

  `flue_loop::network::performance::file` also admits explicit branched tube
  graphs, section-local material gas records, addressed viscothermal sections,
  passive endpoint/interior R-L-C loads, and multiple coplanar radiating
  openings. The jet and windway bind the actual inlet section's complete gas,
  independent of source order. The same physical graph, fractional/loss owner
  and multirate receiver retain the original node addresses through loss
  subdivision. CLI provenance describes every supplied node, section gas and
  selected opening. The single-tube grammar retains its original construction.
  `music_render_wind::flue::network` supplies G0/G3 source/direct-construction,
  gas/branch/receiver sensitivity, callback/PCM and refusal regressions. This
  adds data access to existing physics; it does not add finger-pad mechanics,
  junction end corrections, mutual radiation feedback, gas evolution or
  measured-instrument authority. See `examples/FLUE_NETWORK_PERFORMANCES.md`.

  `FlueOutletPerformance::with_mouth_reaction` selects compact self reaction of
  the actual rectangular mouth in the same exterior gas and area basis as its
  pressure receiver. Optional `mouth_aperture rectangle WIDTH_M LENGTH_M`
  before `mouth_radiation` supplies acoustic dimensions independently of the
  jet's width and cut-up; omission retains the original rectangular window.
  Network observation and compact reaction share this declared aperture.
  Independent aperture dimensions are not yet supported by the DG route.
  This choice does not alter the aerodynamic jet geometry,
  source coefficient or independently selected separation law. The exact
  rectangular static Rayleigh integral and compact
  monopole resistance determine one positive-real `R*s/(s+p)` branch in the
  existing `fs-vfit` midpoint owner. Current-flow pressure and retained history
  enter the original same-sample jet/mouth/duct solve. Original windway, jet
  feedback, gas sources and separately declared contraction loss remain. The
  branch replaces the legacy end correction and authored linear resistance;
  near-field storage and escaped radiation energy have separate ledger terms.
  The file selection is `mouth_reaction compact-baffled` immediately after
  `mouth_radiation`, with zero `mouth_loss_rel`.

  Construction bounds cell work, requires matched static pressure and a compact
  band, checks static quadrature against the exact rectangle, and checks both
  analogue and realized complex impedance and resistance within 5% at 32
  frequencies. These sampled checks are not a continuum or broadband
  certificate. Mouth-to-terminal mutual reaction, free-field body scattering,
  resolved turbulence and measured flute qualification remain outside this
  selection. Physical storage/work and source-to-PCM regressions are authored;
  syntax checks pass, while native compilation and execution remain required.

- `wind_line` module (bead 3ez8g.6.3): the wind articulation runtime —
  `WindLineBank` holds per-fingering exact-FIR characteristic lines
  over a typed `FingeringTable` with the carry-history + crossfade
  switch (the MM bank's D17 lift, plane-wave edition), and the
  char -> VFIT hop: a settled note may leave the exact-FIR phrase
  image for the cheaper rational hold image via the REPLAY-PRIME lift
  (the incoming IIR line is fed the recent outgoing history, outputs
  discarded, before the crossfade), hopping back on gesture
  resumption. Hop-readiness is GATED: entering an image whose
  registry gate is not green refuses structurally
  (`ImageNotGated`). "When to hop" is MEASURED DATA
  (`data/claims/wind-hop-policy.tsv`: settle-detector parameters +
  click-vs-hop-timing table, minted by the fixture, gated on every
  run), not doctrine. Fundamental estimation in the fixtures is
  AUTOCORRELATION — three estimator traps were executed en route:
  a global FFT peak returns the dominant harmonic; even-factor HPS is
  structurally blind on odd-harmonic (closed-pipe) spectra (the
  |X(2f)| factor lands on the absent even harmonic); odd-factor
  products cannot disambiguate f0 from 3f0 (the odd series is
  self-similar under 3x). The period domain has neither ambiguity.

- `piano_vertical` module (bead 3ez8g.5.1): the piano vertical
  composition — still no `Piano` type. Three unison strings (exact-ZOH
  modal images whose mode series carries the stiff-string
  inharmonicity law `f_n = n f0 sqrt(1 + B n^2)`; tension is state) +
  a duplex segment + a small modal board coupled through the
  POWER-CONJUGATE bridge (`F = sum c_k v_k` in, `g_k = -c_k v_b`
  back) + the 87zbd felt hammer island + pedals as coupling states
  (dampers = viscous drag through the strings' own force ports; una
  corda = struck-string count). The Weinreich aftersound EMERGES:
  in-phase unison motion pumps the lossy board (fast early decay);
  the surviving quasi-antisymmetric configuration barely couples
  (slow aftersound). MEASUREMENT LAWS learned executing the gates
  (each an executed misread first): the board-velocity spectrum
  carries the BOARD's own lines — the string law must be read from
  the string state; the board envelope is structurally BLIND to the
  antisymmetric survivors — the two-stage fit runs on the string
  ensemble energy; strong detune SUPPRESSES energy exchange (real
  physics) — the beat gate reads the AUDIO envelope, whose chorus
  modulation persists; and envelope smoothers must be scaled between
  the coupling slosh (~33 Hz, measured) and the detune beat (2-4 Hz).

- vibration gates (bead 3ez8g.7.1): the string/plate/nonlinear
  registry review. The TRUNCATION GATE makes modal honesty executable:
  a row's retained series must have converged (top-retained-mode energy
  share below the authored gate), and the deliberately under-truncated
  falsifier (N=3 at 220 Hz, whose top mode sits mid-band where the felt
  hammer still delivers ~1e-1 of the energy) FAILS the same gate the
  disclosed N=12 fixture passes at 6e-20. The linear<->nonlinear
  SELECTOR THRESHOLD is committed data with fixture provenance
  (`data/claims/vibration-selector-thresholds.tsv`), re-derived every
  run from the KC glide law the fs-nlmodal battery pins at 1e-12 — a
  hand-edited threshold fails. THE GONZALEZ CAVEAT travels with the
  KC/von-Karman promotion: conservation under discrete-gradient
  stepping is structurally blind to gradient errors, so nonlinear
  evidence cites FD-gradient/trajectory oracles, never energy
  conservation.

- Dispersive-waveguide bake-off (bead 3ez8g.7.2): the [F] string
  alternative measured against modal ZOH on ONE shared music-wire card
  with the analytic stiff law `f_n = n f0 sqrt(1 + B n^2)`,
  `B = pi^2 E I/(T L^2)` derived from the card, as the oracle. The
  executed receipt decided KeepBoth: modal keeps every certified claim
  (D21 asserted against the live registry in the gate test); the
  waveguide earns the budget-constrained hero-string scope
  (B_hat/B = 1.014, worst partial 9.8 cents with the linear
  fractional-interpolation phase as the named residual owner, ~10x
  fewer flops per sample, state cost flat in partial count — while
  modal holds FEWER states on this card, recorded honestly). The
  DISPERSIONLESS control on the same loop misses the law by 41 cents:
  the allpass stage, not the loop, carries the stiffness.

- Xylophone bar (bead 3ez8g.12.1): the cheapest FULL filling — a pure
  composition of live parts (analytic free-free bar modes + exact-ZOH
  runtime + a Hertz mallet island; no new physics, only binding and
  gates). THE BOUNDARY CONDITION IS LOAD-BEARING: the free-free
  `beta L` values are the `cosh cos = 1` roots computed in-fixture by
  deterministic Newton (self-verified at machine zero, never
  transcription-trust), giving the non-harmonic ladder
  f2/f1 = 2.7565, f3/f1 = 5.4039; the executed falsifier proves the
  gate discriminates (the PINNED n^2 family measures f2/f1 = 4.000 and
  fails). Material: the Indian-rosewood matdb pack (E_L = 13.53 GPa;
  density = SG*1000 as an authored Estimate, disclosed). The named
  successor is the ARCHED bar (ratio tuning toward 3-4, a
  varying-section solver) — v1 claims the uniform prismatic bar
  exactly and only that.

- Cane reed card (bead 3ez8g.6.2, completing T-Wind): "shaped like a
  reed, made of cane" — a tapered vamp mesh + the CC-BY cane packs
  through the fs-solid reduce lab (f1 = 2709 Hz, the measured
  clarinet-reed class) driving the EXISTING realize_reed_bore loop;
  the authored-parameter path remains the Estimate debug image. The
  LAY is a measured-geometry chart (tip opening + facing length ->
  circular-arc radius) consumed as the rest gap, rolling-stiffness
  hardening the named successor; contact stays the fs-dcontact
  Obstacle pattern. Provenance chain in the card: mesh digest,
  material resolution identity (heel E_L dry; the blank-zone density
  mismatch and the 10x E_T anisotropy DISCLOSED; eta = delta/pi from
  the damping pack), reduction options, retained-compliance
  disclosure. Measured sensitivities: more-open lay +46% threshold;
  harder cane brighter at matched breath (+1.9%) with the ordering
  inverting near the soft threshold (both regimes recorded); the
  stiffness threshold shift is below the ladder rung (logged, never
  overclaimed).

- `glottis` module (bead 3ez8g.8.2): the glottis as the SAME object
  as a lip reed and relief valve — two islands on the menu (D21): the
  1-DOF image with the TITZE SURFACE-WAVE mechanism (the flow sees the
  gap one mucosal propagation time late; a bare quasi-steady one-mass
  valve measured NO phonation to 6 kPa — the delay IS the mechanism;
  wave speed authored Estimate) and the two-mass image (4-state
  coupled pHS with the coupling INSIDE the Hamiltonian, corrected
  Steinecke-Herzel driving, asymmetric seed escaping the unstable
  symmetric rest, dcontact collision). FIVE energy/mechanism defects
  were EXECUTED en route and their causes are recorded in-code where
  each was fixed: the short low-loss tract's FIR wraps non-passively
  (4x-oversampled enforcement, the MM-bank lesson); the junction
  Picard iteration does not contract when Zc U rivals the drive
  (closed-form quadratic solve, the brass lesson); a stale one-sample
  line push adds junction phase slip; the S-H recovery factor's SIGN
  is load-bearing (inverted, the cycle locks antiphase onto the tract
  resonance); and the symmetric rest state is a driving-pressure
  deadlock (S-H's own displacement seed). Fold parameters:
  masses/stiffnesses/length RECEIPTED from the licensed two-mass pack;
  damping, rest gap, collision, mucosal speed authored Estimate — the
  honesty ladder's permitted rung, card-driven either way.

- Glottal flow observation (`glottis::glottal_qois` and the checked
  `glottal_qois_with_options`) uses fixed-window FFT correlation and the
  YIN cumulative-mean normalized difference, with raw-difference parabolic
  lag interpolation. Fractional mean crossings refine coherent periods and
  avoid integer-sample jitter. Pitch-band edges bracket the continuous
  trough; an out-of-band estimate is not clamped into the band or promoted
  to its double period. Constant or unvoiced records identify no pitch
  (`f0_hz = 0`); invalid input is a typed refusal, or NaNs through the
  original convenience API. The pitch band, periodicity threshold, and
  conservative numeric workspace limit are explicit options. DC and up to
  six strictly sub-Nyquist harmonics are jointly fitted by `fs-fft`'s shared
  QR owner. Peak amplitudes and unexplained RMS retain input flow units;
  absent or residual-buried harmonics do not acquire invented finite dB
  levels. Undefined slope or ambiguous mean-crossing jitter remains NaN.
  This is a descriptive, approximately stationary record observer, not a
  clinical voice standard, calibrated voicing probability, or physical
  uncertainty certificate. G0/G1/G3 regressions cover the FFT/direct-sum
  identity, constant/noise/ramp rejection, fractional-period spectral tilt,
  missing fundamental, band edges, Nyquist, physical scaling, and refusals.

- Glottal harmonic observation now refines the YIN-admitted fundamental by
  the shared `fs-fft` bounded local variable-projection solve. The local
  interval is the smaller of 2% of the coarse pitch and one inverse record
  duration on each side, clipped to the declared pitch band. Numerical
  frequency tolerance and the candidate-fit limit are explicit options;
  nonconvergence, nonidentifiability, or resource refusal propagates through
  the checked API instead of publishing an unfinished estimate. The observer
  reuses the refined joint fit. If the refined pitch admits another member of
  harmonics 1..6 below Nyquist, a second refinement includes that member under
  a tightened Nyquist bound and the remaining candidate-fit budget. The
  2% neighborhood admits at most one such extra order. Actual crossing
  jitter remains independent of this stationary frequency refinement. The
  missing-fundamental regression retains positive flow, h2/h3 dominance and a
  10,000:1 amplitude range down to h6=1e-8 in flow units; an omitted
  inharmonic tone remains in the physical residual. This local correction is
  not a global pitch search, octave certification, or uncertainty estimate.

- `tract` module (bead 3ez8g.8.1): the tract as A(x,t) with LICENSED
  charts (Assaneo 2011 PLOS CC-BY Table 2 — the authors' own 10-tube
  /a/ and /u/; Spanish, acoustically-inverted, formant-valid register
  disclosed; Story 1996 is unlicensable and the Dresden CC0 meshes are
  the named MRI-grade upgrade vein). D22 as a WHOLE-CURVE gate: char
  impulse-spectrum envelope vs the TMM |Z_in| envelope of the SAME
  walled chart. STRUCTURAL FINDS the lane earned: characteristic_line
  consumers must pass VOLUME-normalized zc (rho c/S — an area-free
  rho*c flattens R to +1 and silently erases the tract shape; the
  glottis module carried this defect and its gates were re-based under
  the corrected load); a rigid-wall short tract realizes ACTIVELY (its
  ringing outlasts the FIR window and time-aliasing makes |R| > 1 —
  the tissue walls are required for the REALIZATION, not just the
  physics); a uniform passivity clip that fires flattens the shape
  (backstop only); dense realization via characteristic_line_dense
  (n_fft override). The morph carry needs p-minus CONTINUITY or every
  rebuild micro-clicks.

- Vowel gates (bead 3ez8g.8.3, completing T-Voice): the first voice
  claims — the vowel from the licensed area function, the source from
  the self-oscillating folds, their INTERACTION through the load at
  the glottis. The corpus loop closes (rows
  acoustic-assaneo-2011-tract-{a,u}): the charts' TMM F1 lands in the
  published Spanish vowel CLASS bands (827 / 340 Hz; a class check,
  never a per-subject match) and the length-shifted falsifier fails
  the same gate at 1104 Hz. THE FLAGSHIP EMERGENT GATE: sweeping the
  fold frequency toward F1(/u/), the phonated f0 rides ~35% above the
  free-fold frequency (load stiffening) and REGISTER-BREAKS at the F1
  ceiling (deviation -147 Hz at the top of the sweep) — the physics a
  source-filter vocoder structurally cannot have, measured. The
  1-DOF-vs-two-mass bake-off receipt is EXECUTED (KeepBoth:
  cheap-ensemble vs solo-realism; slope -6.6 vs -11.5 dB/oct on the
  identical fixture); the ooh/aah ladder is minted for the owner's
  ear.

- `pickup` module (bead 3ez8g.9.2): the Faraday pickup as a
  PURE-READER PORT — D20 as structure (pose + per-mode window gains
  only; `emf()` reads the string owner's modal velocities; no step, no
  state, no re-discretization anywhere). The raised-cosine aperture
  window is unit-normalized so the point-pickup limit is EXACT; the
  window is what makes bridge-vs-neck voicing physical. v1 provenance:
  linearized flux with an authored dipole-like height profile
  (Estimate, disclosed); the optional material field bake upgrades this
  provenance without changing the runtime's Faraday operator.
  `pickup::field` supplies that opt-in field-backed path: original planar cell
  material states, fixed magnet orientation, signed winding turns and declared
  moving mesh vertices feed `fs-em`'s existing fixed-remanence winding solves.
  The bounded Hermite linkage map uses h/h2 derivative comparisons and fresh
  midpoint field checks with explicit SI tolerances and field-evaluation limits.
  All field samples retain native residual reports; these sampled comparisons
  are numerical estimates, not continuum certificates. Invalid deformation,
  moving winding/magnet cells, unsupported gaps and cancellation refuse.
  `MaterialStringPickup` reads the existing homogeneous material-string owner's
  mass-normalized state without allocation, field solves or a second integrator.
  Sensed cells must share its material card and query state and translate rigidly;
  material/geometry or basis replacement invalidates the binding, while gas-only
  updates remain visible through the owner's motion. The generic map retains
  the original Faraday observer and its reciprocal discrete-work force port.
  This path is open-circuit magnetostatics: no static magnetic attraction,
  electrical loading, induced-current losses or saturation is inferred. Planar
  geometry and aperture reduction still require physical validation.
  G1/G3/G4 `tests/pickup_field.rs` checks an independent exact layered field,
  remanence/permeability data substitutions, reciprocal work, native-string
  observation, stale-source refusal, numerical budgets and cancellation.

- `speaker` module (bead 3ez8g.9.4): the Thiele-Small driver + sealed
  cabinet as ONE port-Hamiltonian — the Bl gyrator is an antisymmetric
  J pair (power-exact by construction), the sealed box is
  `Sd^2 rho c^2 / V` extra cone stiffness (exactly why the classic
  resonance-shift algebra holds), and one fs-plate-class panel mode
  couples through the shared cavity term (the cabinet sound as
  physics; X-Consist with ONE retained mode, disclosed). Excursion
  past the card's x_max REFUSES — the suspension nonlinearity is not
  modeled and silent linear forte is the lie the refusal prevents.
  Card = authored datasheet-class Estimate; low-ka piston load
  disclosed; cabinet IR packs are FORBIDDEN and absent.
  `speaker::harmonic::BaffledTsDriver` is the continuous-frequency companion:
  the shared full Bessel/Struve piston impedance loads the original coil,
  suspension, optional cavity and panel equations at each actual frequency.
  Its reactive part supplies the fluid inertia once; neither the transient
  owner's fixed radiation resistance nor its compact added mass is retained.
  The original `fs-phs::harmonic` solver owns resource admission, cancellation,
  equilibration and equation residuals. Peak phasors use exp(-i omega t).
  Electrical supply, coil/suspension/panel heat and exported acoustic power
  are independently reported; pointwise equivalent inertia is not published
  as stored acoustic energy or promoted to a causal transient load.
  Cone pressure uses signed circular-piston directivity, half-space normalization
  and physical propagation phase. Caller limits bound the far-field kernel's
  range/phase approximation; no relative-pressure bound at a null is claimed.
  G1/G3 tests compare an independent circuit reduction, low/high-frequency
  piston limits, panel antiresonance, reciprocal power, and hemisphere-integrated
  acoustic intensity. This bounded linear infinite-baffle model supplies no
  panel radiation, cone breakup, mutual exterior loading, room, nonlinear
  excursion or experimentally validated driver data.
  `TsDriver::new_with_media` and `MagneticTsDriver::new_with_media` retain
  distinct exterior and sealed-cavity gas snapshots. Exterior gas owns the
  original radiation load and pressure observer; cavity `rho*c^2/V` owns
  compression and reciprocal cone/panel coupling. A box needs an explicit
  cavity gas, and both static pressures must match: pressure-bias equilibrium
  is not modeled. The single-gas constructors retain identical-medium behavior.
  `ThermalMagneticTsDriver::from_material_with_media` uses that same assembly;
  coil-heating trials preserve both gas states without inventing gas heat flow.
  The opt-in `MagneticTsDriver` retains that original mechanical assembly and
  admits caller-supplied quadratic Bl(x)/L(x) profiles over an explicit closed
  displacement interval. L must remain positive throughout, including interior
  extrema; the interval must fit inside the suspension's original travel limit.
  Total linkage lambda is stored with H_m=(lambda-psi_pm(x))^2/(2L(x)), where
  psi_pm is the analytic integral of Bl. The same gradient supplies current,
  Bl(x)i + L'(x)i^2/2 cone force, and reciprocal voltage including L'(x)i*x_dot.
  No additional linear gyrator duplicates the interaction. The existing
  descriptor discrete-gradient solver owns the common step and work/heat ledger;
  rejected trials do not change accepted state. This nonlinear path is second
  order and needs timestep refinement; the linear TsDriver keeps ExactZoh.
  Tests in `speaker/magnetic/tests.rs` cover constitutive energy derivatives,
  independent voltage/power identities, finite-domain admission, constant-model
  convergence, nonlinear cabinet/panel evolution, and source/heat balance.
  These bounded profiles do not add data authority, current-dependent core
  saturation, hysteresis, eddy losses, coil heating or nonlinear suspension.

- Electric chain (bead 3ez8g.9.5): `tests/electric_chain.rs`
  composes string -> pickup -> tone RC -> triode island -> TS driver
  with every stage a named, receipted object. The composition laws:
  ONE string owner (a bare string stepped alongside stays BITWISE
  identical — the chain reads, never owns); voicing survives the
  chain (per-harmonic bridge/neck output ratios = pickup window-gain
  ratios on string-carried lines); the harmonic ladder is
  DEVICE-made (single-mode string, so 2f0/3f0 at the speaker is
  constitutive curvature); the box composes transparently
  (chain-measured excursion-per-volt = driver-alone). The
  split-vs-full-DAE bake-off receipt
  (`tests/receipts/electric-split-vs-full-dae.bakeoff`, KeepBoth)
  pins the split image to the analytic gain and to the monolithic
  1e-14 authority. MEASUREMENT LAW recorded: harmonic line
  projections need INTEGER-CYCLE windows — at a non-integer cycle
  count the fundamental leaks ~1.7% into every harmonic line,
  drive-independent, and read as a 27% violation of a ratio law
  that holds.

- Composed guitar body (bead 3ez8g.7.4): `tests/guitar_body_gate.rs`
  gates the braced-top + ladder-braced-back + Helmholtz model's
  BRIDGE-ADMITTANCE triad against the Carcagno 2018 BR-guitar corpus
  rows (F1 +5.9%, F2 +9.9% in 15% envelopes; F3 -32.6% in an
  authored 40% envelope with the mechanism DISCLOSED — the flat
  rectangle surrogate has no X-bracing/arching, so its third
  breathing family sits low: the named weakest claim). Executed
  falsifiers: the midline-brace order flip (frequency-order breathing
  selection picks a 6e-15-volume mode; select by max |INT phi dA|)
  and the MEASURED refutation of nodal-line folklore (an ECCENTRIC
  brace on the nodal line shifts the mode 11.6% via torsion +
  eccentric-membrane coupling). The wolf note is a pHS pair: tune
  EFFECTIVE frequencies (k = m w^2 - k_c; the coupling spring detunes
  a naive pair by 10% and no hybridization happens), then the string
  on the breather beats at DETRENDED depth 0.046 (the estimator
  lesson: divide out decay before measuring modulation) vs 0.990
  clean off-resonance. Truncation sensitivity is PER-FREQUENCY
  (spread 23.6x) and its direction is fixture-specific: here
  between-peaks suffers more than on-peak because low-index resonant
  modes survive basis halving.

## Conformance tests

`src/broadband_radiation.rs`: G0 complex/real-SH round trip/reconstruction against
`fs-bem`; admission, neutral replay, properness/direct-feedthrough retention, and
held-out mutation; causal transactional runtime superposition and bitwise replay.

`src/vibroacoustic.rs` unit tests: closed-form overlaps, two-oscillator
split + dropped-coupling mutation, added-mass first-order falsifier,
exact energy balance + one-sided normal-flip mutation, reciprocity,
truncation refusal, cavity/resonator pins (incl. the elongated-cavity
regression), multi-mode pencil vs independent determinant, refusals,
determinism. `tests/vibroacoustic_casebook.rs`: the fs-plate ->
fs-couple -> fs-bem composition — box-with-flexible-top vs independent
perturbation with the sealed-cavity stiffening direction, coupled FRFs
with observed per-frequency truncation deltas, and the radiation solve
whose per-period energies audit GREEN through a real
`WindowAuditReport` (with the `entropy_on_nonthermal_power_port`
refusal retained as doctrine: dissipation leaves through THERMAL
ports); plus the guitar T(1,1) comparison against Carcagno et al. 2018
with the exact two-DOF product invariant and the material-vs-measured
damping inequality.

`tests/couple.rs` (38 cases): v2 scalar-seed migration goldens for all three
legacy kinds and bitwise whole-result migration of the added-mass fixture;
rotational/electrical/magnetic/chemical watt-dimensional goldens, required-role
admission, raw-oracle refusal, and semantic kind mismatch refusal; schema
fail-closed metadata; localized junction mismatch; non-scalar refusal by the
scalar evaluator; field-duality measure dimensions, density-side assignment,
measure-application overflow, and distinct effort/flow spaces; all four
primitive descriptors and identity-alias refusals;
complete moving-enthalpy stream admission and canonical constituent ordering;
one proved explicit-chemical contribution;
bit-exact enthalpy/internal-energy pressure-work equivalence; one-ULP and
foreign-context crosswalk/stress-evidence refusal; wrong-identity and
non-positive-density refusal; all four exact Euler/Legendre coordinates
reconstructing the same canonical enthalpy and chemical double-count refusal;
proof-bound axis mismatch and species-mode/element-axis refusal;
empty/duplicate/non-finite/wrong-orientation stream refusals and
declared-energy mismatch;
green first-law/mass/element/charge/entropy closure with explicit signed
advected/diffusive/heat/radiation rows and `T0*S_gen`; retained red first-law
and negative-entropy reports; exact missing-row and foreign-evidence-role
refusal; exact manifest/element/charge evidence-subject refusal; same-interval
energy/accounting-context rebinding refusal; stream/chemical element-schema
and typed power-axis refusal, including cancelling forbidden entropy slots,
inexact constant-temperature energy, and constant-temperature radiation;
boundary-temperature and nonfinite fail-closed admission; all six
permutations of a cancellation-sensitive compensated sum; an explicitly closed
isolated entropy window; one-sided internal-pair refusal; and a hidden internal
entropy source localized to the exact two rows, ports, and exchange while an
external row masks its pair-local energy defect from the global first law;
legacy power-conjugate ports; exact interface power conservation and
incompatible-port refusal; energy-audit imbalance and non-finite alarms; the
Aitken Δ² factor; naive staggering diverges where Aitken stays stable; Aitken
accelerates over stable fixed relaxation; light added mass converges naively;
determinism.

`src/mm_line.rs` `mm_line_tests`, cases ml-001..ml-006 — held-fingering
line fidelity vs the MM authority (0.00 cents worst over three valve
combinations; the scattering LOOP closes the reflection at the
junction's own sample and reads 2.86 cents worst, gated at 4), the
realization bake-off receipt, the valve
no-click gate (carry+crossfade within 3x the control's slew; the
hard-cut falsifier is at least 2x worse), per-line passivity on the
stored taps, the cup shunt vs its analytic transform, and
causality/splice/refusal arms (the causality gate is the executable
witness of the e^{-i omega t} conjugation discipline).

`src/brass_loop.rs` `brass_gate_tests` + `brass_observer_tests`, cases
bl-001..bl-008 — the EMERGENT gates: lock 238 cents above the
cup-loaded column peak (outward-striking band) and 281 cents from the
lip alone at 10.4 kPa; tension walks the lock 208.9 -> 445.1 Hz across
distinct column peaks; valves drop the whole series (-125/-200 cents,
monotone, no pitch anywhere in the inputs); slotting (+78 cents of
continuous lipping inside a slot, 1310-cent jumps between); the cold
horn plays -50.4 cents flat vs the -48.6 sound-speed prediction with
temperature entering ONLY through the gas card; the no-click valve
change (envelope held through the switch; the fresh-voice hard-cut
falsifier collapses to 40% RMS); structural honesty + refusal arms;
and the radiated-observer demonstration from the SAME bake's SH
directivity table at the lock frequency (compact-mouth near-omni is
the correct physics at ka 0.09).

`src/wind_line.rs` `wind_articulation_tests` + `wind_hop_tests`,
cases wa-001..wa-004 + the committed-policy gate — the sigma slur
(+265 cents with the envelope held through the change), the register
vent (MEASURED, not the textbook sketch: with an INERTIVE 1.8 mm /
6 mm-chimney vent the twelfth-class flip happens at 0.70L, x3.36,
18 cents from the TMM authority's own peak, while the 0.30L vent
leaves the lock in the low regime — position matters AND the lock
follows the impedance peaks), the structural gate refusal, and the
hop policy (settled hop click 0.504 vs the early-hop falsifier 0.864;
the note survives on the hold image; the committed policy artifact
re-gated every run with the measured nuance that the settle detector
is conservatively late — the right side to err).

`src/piano_vertical.rs` `piano_vertical_tests`, cases pv-001..pv-005 —
THE TILT CONTRAST on the composed vertical (felt trend 159.7 vs the
matched linear spring's 1.000: the hysteresis is audible, not
asserted); the inharmonicity law measured in the audio (worst 4.23
cents over partials {1..5,7} with the B=0 falsifier missing by 15.04;
partial 6 disclosed as the duplex collision line, partial 8 as
weakly-excited at the strike point); the two-stage decay (string
ensemble 3.28/s early vs 1.15/s late) with audio-envelope beats
scaling with detune (0.370 s -> 0.239 s at doubled detune); pedal
topology (dampers collapse the note; una corda energy ratio 0.67 with
the third string still fed through the bridge); and the dissipation
ledger + bitwise replay (worst per-window growth 0.0).

`tests/vibration_gates.rs`, cases vg-001..vg-003 — the truncation gate
with its executed under-truncation falsifier (6.2e-20 vs 1.0e-1 across
a 2e-3 gate), the selector-threshold re-derivation pin (relative 1e-8
against the committed artifact), and the gate-summary enumeration
(every vibration-filling registry row must appear in
`data/claims/vibration-gate-summary.tsv` with a status MATCHING the
registry; the plate 1-port's honestly-ungated row is asserted by
name).

`tests/bakeoff_dispersion.rs` — the committed-receipt gate (bytes
decode, hash, KeepBoth outcome, D21 modal-row-untouched against the
live registry, the newcomer's own row present) and the
passivity/control gate (loop bounded and decaying under g < 1 with the
unit-magnitude stage; constructor refusals; the dispersionless
control's partial 8 sits within 20% of harmonic while the law demands
+18 cents).

`tests/xylophone_bar.rs`, cases xb-001..xb-004 — root self-verification
(residual at machine zero, classical table at 1e-8, asymptote check) +
named refusals; the rendered ratio gate within 1% with the executed
pinned-family falsifier; strike-position physics emergent (mode-1 node
found in-fixture at 0.2242 = the cord-mount point; striking it kills
the fundamental ~8 orders, the center strike kills partial 2, and the
non-noded partials keep speaking); and the listening chain gate
(WAV -> sidecar -> receipt, Unadjudicated).

`tests/cane_reed_card.rs`, cases cr-001..cr-003 — the provenance-chain
card mint (identity reproducibility, retained-compliance floor, the
reed-class f1 band, the lay's facing-scale radius); the sensitivity
sweeps with a SUSTAINED-oscillation speak detector (tail at 2% of
blowing pressure and non-decaying — the absolute-threshold detector
fired on ringing at every rung, twice, before this); and refusals by
name (a card without a modulus cannot resolve; lay-chart nonsense;
an empty orifice face).

`src/glottis.rs` `glottis_tests`, cases gl-001..gl-005 — both islands
self-oscillate against the live tract with onset curves in the
physiological range (1-DOF 188-276 Pa, two-mass 404-865 Pa, rising
with stiffness, none on the ladder floor — the detector demands
DETRENDED sustained output after a settling ramp read as phonation on
the first run); the collision/damping ledger (undriven kick decays,
per-step dissipation non-negative); the two-mass vertical phase
(upper lags lower 0.306 rad at the fundamental) with faster-closing
flow; the bake-off FIXTURE QoIs on the identical tract + card (f0
152/173 Hz, OQ 0.78/0.82, slope -6.6/-10.0 dB/oct — the receipt is
vowel-gates' to mint); refusals by name.

`src/tract.rs` `tract_tests`, cases vt-001..vt-004 — the whole-curve
D22 gate (envelope correlation 0.87//a/ and 0.90//u/ with dominant
resonances within 50 Hz; identity /a/ >> /u/); the yielding-wall F1
raise (340 -> 401 Hz, nearest-peak comparison past the wall's own
~100 Hz mode); the u->a morph (100 Hz control rate, carried, clicks
18x less than the five-cut waypoint control on the REFLECTED wave,
F1 track 353 -> 806 Hz logged); licensing-first refusals (unlicensed
data refuses by name). The glottis module's gates were RE-MEASURED
under the corrected load in the same change (two-mass onsets
404-715 Pa physiological and monotone; the 1-DOF's ~5.8 kPa
non-monotone onsets recorded as the fixed-delay mechanism's real
behavior; flow skew logged, not asserted).

`tests/vowel_gates.rs`, cases vw-001..vw-004 — the formant class gate
with its executed length-shift falsifier (a uniform AREA scale would
not move formants; the falsifier must shift LENGTH); the sung
interaction sweep (ten stiffness points crossing F1 with per-point
f_free/f0/deviation logging); the committed bake-off receipt gate
(KeepBoth scopes + both rows alive in the registry); and the
vowel-ladder listening chain (Unadjudicated).

`src/pickup.rs` `pickup_tests`, cases pk-001..pk-005 — station
voicing with every per-mode EMF line matching an INDEPENDENT
Simpson-integrated analytic gain within 5% (bridge centroid 452 Hz vs
neck 328 Hz); two mixed pickups combing per the gain-sum prediction
with an executed notch; the exact point-pickup limit (5e-4) plus the
pure-reader law (identical repeated reads; EMF recomputable from the
public gains and the same state slice); refusals by name; and the EMF
driving the .9.1 circuit DAE end-to-end (RC-load rolloff matching the
first-order analytic response at 0.25x and 4x cutoff with supply
defects at 1e-21).

`src/speaker.rs` `speaker_tests`, cases ts-001..ts-005 — the sealed
box shifts resonance per fc = fs sqrt(1+Vas/Vb) (measured 64.7 vs
65.1 Hz analytic); the electrical impedance matches the analytic TS
motional curve to the printed digit at four frequencies; the panel
colors the response near its resonance and NOT far away (executed
contrast against the rigid box); a 60 V resonance drive trips the
excursion refusal while a gentle drive never does; refusals by name
plus the supply-defect ledger at 3e-19.

`tests/electric_chain.rs` (7 + 2 mints): ec-001 compose + logging +
the one-owner witness; ec-002 voicing ratios within 6% (gain-floored
to string-carried lines); ec-003 H2 0.07% -> 0.86% across grid drive
0.03 -> 0.40 V RMS with H2 >> H3; ec-004 composed-vs-direct box
transfer within 0.04% at 65/260 Hz; ec-006 committed bake-off
receipt re-asserted (split gain 59.950 = authority 59.950, analytic
59.923; pluck L2 4.6e-12); ec-008 listening digest chain
(electric-clean-vs-driven, Unadjudicated); ec-009 the electric
gate-summary artifact enumerates every registry row.

`tests/guitar_body_gate.rs` (5 + 1 mint): gt-001 triad + product
invariant (defect 0.0) + damping honesty (material-only Q 39/54/78
above measured 34/18/36 as it must be); gt-002 the executed
midline-brace negative control; gt-003 wolf note; gt-004 truncation
spread; gt-006 guitar-strum listening digest chain (Unadjudicated).

