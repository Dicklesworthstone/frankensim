## Nonlinear membrane damping and source-bound finite heat storage

`membrane_thermal::advance_candidate` composes the existing nonlinear
`DampedMembraneMotion` with `DepositionSimulation::advance_observed`. One
unique active caloric receiver belongs to each original reference triangle.
Its `SheetCaloric` retains the identical complete original sheet source and
must be the exact caloric law already borrowed by that thermal cell. Volume is
original area times reference thickness; mass is the original P1 corner-density
integral, including nonuniform reference inventories and constrained coordinates.
An initial 512-epsilon relative allowance covers arithmetic only. Subsequent
cell/material/body/segment, geometry, mass and caloric-instance checks are exact.

The native thermal owner chooses every conductive/event/source interval, also
bounded by the mechanical step ceiling. Analytically integrated per-triangle
Kelvin--Voigt heat becomes average power on that same interval. Original source,
conductance, ambient/reservoir, pairwise-transfer, transition and observer hooks
remain active. No instantaneous-Cp temperature increment, second thermal solver,
numerical-defect heat or support-impulse heat is introduced. The source chart
performs the original nonlinear enthalpy-to-temperature conversion.

Mechanical state is borrowed until the first actual interval, then cloned once;
its prepared inertia and damping buffers are reused across substeps. The thermal
owner supplies whole-call transactional cloning. Motion, thermal inventory,
observer history and energy records publish together after final checks;
refusals and cancellation expose none of their half-steps. Compensated accounts
retain physical heat, held-load work, support impulses and signed numerical
defects. Absolute substep mechanical defects and per-cell heat-handoff errors
have whole-call budgets, preventing cancellation of opposite errors.

The combined physical residual is mechanical storage change plus native thermal
storage change minus incoming inventory and external work/heat. It remains an
uncorrected numerical residual; native mechanical defect is separately disclosed,
never deposited. Frozen elastic/viscous coefficients, natural reference and
geometry must stay within the caller's explicit fully solid temperature interval.
This is one-way finite self-heating, not evolving-temperature material feedback,
reversible thermoelastic work, phase-dependent mechanics or a stability/accuracy
certificate. Native transport order and nonlinear split error require refinement.



## Explicit rigid-impact loss heating in original finite caloric cells

`friction_heat::directional::impact::heat::apply_source_rigid_elliptic_impacts_with_heat`
composes the original simultaneous directional impact owner with
`DepositionSimulation::apply_heat_j`. A unique-cell-ID map explicitly assigns
nonnegative fractions of the one GLOBAL native impact dissipation. Fractions
must total at most one; compensated admission retains even sub-ulp overspend.
The remainder stays in `ImpactHeatReceipt::unallocated_energy_j`. This policy
does not infer a heat-conversion fraction or attribute total collision loss to
individual contacts, normal/restitution loss, friction, wear or radiation.

The original impact/source solve executes once using incoming conditions and
the complete native body/contact mobility. Its dissipated joules reach the
existing finite cells through their actual masses and original nonlinear
caloric laws. Cp/H, Cv/U and latent coordinates retain their source identity and
storage basis. An H receiver retains its pressure-work convention; the combined
account does not relabel H or declared sensible storage as internal energy.
No snapshot heat capacity, impulse-derived pressure, contact area or event
duration is introduced. Empty/zero allocation retains the original mechanical
endpoint and leaves thermal inventory unchanged.

The event occurs at the receiver's accepted physical time; the caller owns
mechanical event alignment and must process due thermal events first. Original
geometry, contact links, source watts and the thermal clock are preserved.
Native inventory limits and slice consistency are checked before staging.
Bodies, source histories and thermal state remain candidates until caloric,
partition, combined-work and final cancellation checks pass. External effects
inside caller callbacks remain the caller's responsibility.

The combined residual sums the original mechanical defect, original caloric
storage defect and explicit partition rounding defect:
`Delta K + Delta E_original + unallocated_loss - support_work`.
Each upstream owner retains its own bounds. No residual is credited as physical
heat, and represented cell-energy errors cannot cancel away inside the original
thermal audit. The partition helper allocates one map per event; no extra
thermal solver, physical time integration or constitutive model is added.

Three native tests execute the exact production partition/deposition helper
against the original thermal crate: distinct analytic Cp(T)/Cv(T) and latent
storage with fixed time/identity, fraction and representability limits, and
late cancellation/domain refusal with deterministic retry. Three full native
mechanical integration regressions are included for an independent 1.36-J
collision, empty/zero mechanical equivalence and joint body/history/caloric
rollback, but were not executed here because external workspace dependencies
are absent. The complete thermal production library compiles; its broader unit
test target separately remains blocked by two pre-existing advection observer
fixtures missing required `event`/`step` methods. No full-workspace verification
claim follows from the focused helper execution.

## Evolving native joint motion with finite-cell friction heat

`friction_heat::joints::trajectory::advance_joint_friction_heat_trajectory` composes the original constrained native half steps and `propose_joint_friction_heat` over a nonempty, strictly increasing, caller-authored time grid. Actual native midpoint states and physical time enter the temperature/material source callback. Every interval recomputes joint geometry and solves friction strength, native recoil, interface temperatures and original finite-cell heat together. Only the native rigid owner advances poses.

The original thermal clock determines each represented interval. Known drive/load switches belong at grid endpoints. Controls bound steps, interval duration, aggregate nonlinear trials, native/thermal inventories, absolute accumulated reaction work, signed total energy residual and the sum of absolute interval energy residuals. Only final body/thermal/source receipts and at most six cumulative joint/link work accounts are retained; memory does not grow with the number of steps. Late failures or cancellation publish none of the intermediate candidates.

Cumulative work accounts distinguish original wrench work, axial drive work, continuous friction loss, converted heat, explicit nonthermal exports, outside heat and measured bilateral residual work. Endpoint mechanical change uses paired native momentum/velocity increments with the original full inertia and gravity potential; caloric change uses each original cell's mass and stored coordinate. No balance residual is added to heat. After a phase with equal actual receiver temperatures, that caloric temperature seeds the next interface solve; it remains a Newton guess, avoiding spurious tiny energy corrections when friction stops.

The thermal routes remain explicitly reduced, fixed-geometry bearing/rail receivers; no moving footprint, material remap, surface recession or flash-temperature law is inferred. Backward-Euler heat makes this a first-order thermal composition despite symmetric mechanical splitting. Caller time refinement remains required: exact static balance belongs in joint drive ports, while separately split body wrenches can introduce artificial axial slip/heat at finite step size. A small energy residual is not a trajectory/heating accuracy certificate.

Four direct regressions cover event-aligned stopping/reversal poses with 1/2 J and 7/18 J loss, eight-step temperature-dependent friction against an independent quadratic recurrence, an authored drive switch with gravity/wrench/outside-heat and nonthermal-export accounts, and late source/cancellation or malformed-grid atomic refusal.

## Hardness-limited sliding, wear and finite caloric feedback

`friction_heat::plastic_wear::PlasticSlidingHeatContact` binds an accepted
`SourcePlasticSlidingWear` checkpoint to its original ordered material cards,
source-resolved Cp/density charts and an existing `FrictionHeatRoute`.
`step_candidate` solves original cell-energy increments and the massless
interface temperature using the existing `FrictionHeatStep`,
`FallibleFiniteDifference` and Newton--Krylov owners. No parallel thermal
balance or nonlinear solver is introduced.

Each residual evaluation, directional probe and final replay starts from the
same accepted source and complete indentation, shear and wear history. The
original caloric objects give separate body temperatures; the normal/wear
temperature factory refreshes original elasticity and each body's pressure
hardness at those temperatures, while original friction and Archard coefficients
follow the solved interface temperature. The selected normal boundary is held
approach, maintained force or a held finite-stiffness device. Native Coulomb
motion uses the resulting actual cell pressures, and only its solved slip
length drives each body's original Archard removal and subsequent normal
re-equilibration. Thermal coefficients are endpoint values; wear remains the
native beginning-pressure split, requiring physical time/motion refinement.

The original finite-cell geometry and endpoint conductivities determine signed
heat into both bodies and a solved common interface temperature. No prescribed
heat fraction or arithmetic temperature mean is inserted. The source callback
owns original conductivity claims and other transport, retaining its query and
selection evidence in the returned payload. The route excludes its old direct
link to avoid duplicate conduction. Only actual native Coulomb dissipation
enters this heat source. Normal plastic work and signed temperature/recession
parameter work remain separate; their explicit external-work budget is the sum
of absolute parameter contributions. Native contact identities and the thermal
cell, interface and aggregate balances have independent acceptance tolerances.

Admission retains the original receiving body IDs, cell geometry, material IDs,
masses, caloric objects and clock. Original caloric and normal sources must share
ordered card identities, complete fixed conditions and quantity axes, and the
selected density claim. Accepted bulk source temperatures must equal the
thermal checkpoint; melted material refuses. The original active grid-cell
mask must have the thermal face's area and remain unchanged after both
temperature refresh and wear. Pressure may redistribute within that mask.

This is an explicit fixed-footprint, fixed-thermal-inventory approximation.
Cumulative removed volume is bounded against each ORIGINAL receiving cell
volume, never a renewed per-step denominator; zero permitted removal is valid.
The fraction bounds an omitted effect, not its thermal error or experimental
validity. Removed material has no modeled mass, enthalpy or debris export.
Reciprocal thermoelastic free energy, conversion of normal plastic work to heat,
moving contact coverage and storage in caller-supplied tangential body response
are not closed by this owner. Its local contact identities and fixed-inventory
heat balance imply no complete body-energy, entropy or passivity theorem.

The candidate seals successor source, contact history and thermal state for
joint adoption. Bounded trials, Newton iterations and thermal inventory,
source-domain or footprint refusal, excessive cumulative removal or parameter
work, nonconvergence and cancellation publish no partial physical history.
The report retains both native mechanical operations, source callback evidence,
friction heat, normal/plastic/parameter accounts and cumulative removed volumes.

## Static breakaway and kinetic sliding through joint heat and motion

Existing joint-friction thermal source callbacks now accept explicitly authored
static/kinetic native port variants. Every nonlinear temperature trial resolves
the selected native Coulomb/stiction law and supplies its continuous friction
work to the original finite thermal network. Holding generates no artificial
friction heat; released and already sliding joints use kinetic capacity. Static
capacity, kinetic capacity and conductivity may be resolved at the same
candidate interface temperature. Ambiguous mechanical continuation remains a
refusal, and accepted inputs remain borrowed until a complete candidate passes.

Three cross-module tests distinguish hold from slide under identical load,
solve temperature-dependent sliding friction against a finite two-cell thermal
response and an independent closed-form balance, and evolve an authored motion
grid through hold, breakaway and continued sliding below the static limit.
The last case checks actual poses, endpoint rates, drive work, kinetic change,
integrated friction loss and thermal storage together. No new capture speed,
implicit presliding history, inferred heat fraction or imposed thermostat is
introduced.


## Scenario pressure histories on native fluid reservoirs

`scenario_pressure` connects declared `IncompressibleFlow/PressureOutlet`
boundaries to the original specialized D3Q19 and generic D3Q19/D3Q27 open
boundary owners. `PressureReference` explicitly distinguishes signed gauge
pressure from positive absolute pressure with a supplied positive reference.
Both use `MaterialScaling::density_lattice_from_pressure_perturbation`; the
adapter does not infer an acoustic equation of state or a pressure-work model.

Uniform values and `TimeSignal` histories retain SI dimensions, the original
interpolation/clamping rules and cancellable validation/evaluation. Scalar
pressure is invariant under a checked frame transform. Spatial profiles require
an explicit geometry binding and refuse this uniform-time lane. Selected native
pressure profiles also refuse, including unchanged dormant uniform targets;
unselected controls and profiles remain active.

A six-face batch prepares every target before one native rebind. Sampling,
reference, units, density, native-kind or cancellation refusal preserves all
controls and fluid state. This updates held controls; native stepping retains
collision, force correction and reservoir transport. G0/G3 tests in
`tests/scenario_pressure.rs` exercise both native quadratures, actual boundary
density after steps, gauge/absolute equivalence, histories, profile precedence,
transactional refusal and cancellation. No new dependency or solver is added.

## Changing contact area with physical conduction and exposed cooling

`FrictionHeatStep::evaluate_with_area` and its entropy variant stage mechanical
area first, then use the native contact-only implicit endpoint ingress before
resolving material transport. The callback sees actual covered/exposed surface
areas and a suppressed direct link. Half-cell fluxes use original conductivity
times the new area divided by the original distance; conductivity itself is
never rescaled to imitate geometry. The fixed and changing-area APIs share
their physical transfer, energy and entropy assembly.

The accepted thermal state and `FrictionHeatAreaTrial` successor route retain
the same new coverage and source-derived series conductance. Actual ambient
cooling follows native coverage without a compensating source. Areas must be
positive and fit both original prismatic faces and aggregate surface coverage.
Complete separation refuses because a vanished massless interface has no
unique temperature. Cells, masses, caloric objects and geometry remain fixed.

This is a uniform-face endpoint discretization requiring time/spatial
refinement, not lateral heat spreading or a swept-contact model. Discrete
active-cell changes may prevent a nonlinear root from converging. Three
focused regressions specify independent shrink/growth energy and entropy
balances, callback-visible coverage and source-k rebinding, exact-target event
order, and transactional late failure including a caught contact-update error.

## Spatial thermal networks in impact pressure playback

`render::plate::impact::audio::network::ThermalNetworkSource` retains an original
heated or liquid-film impact source, its `DepositionSimulation`, thermal observer,
reservoir/cell bindings, step controls and contact-material resolver. Its
`ImpactSource` implementation calls the original `step_with_thermal_network`
once per mechanical sample. It does not first take an adiabatic step, reconstruct
caloric properties, or deposit mechanical dissipation a second time.

`ImpactPressureRenderer` therefore observes the accepted network-coupled
trajectory through its existing radiation and propagation histories. The same
source can pass through causal decimation, pressure mixing and WAV streaming.
`last_mechanical_frame` retains the original `NetworkImpactFrame`, including
thermal substeps, external heat and the combined mechanical/thermal account.
Source, thermal cells and observer are available read-only between callbacks.

The original coupled owner admits complete bindings, matching physical clocks,
caloric identity, source domains, thermal and mechanical budgets, and cancellation.
A rejected physical sample publishes none of those states; the pressure host
retains its existing failed-callback poisoning contract. External effects in a
captured resolver are not rolled back. Thermal splitting, fixed geometry,
unloaded-only contact refresh and one-way acoustic observation keep their
original limitations. Construction retains inputs without advancing or admitting
a physical step. The focused regressions compare direct/network-observed states
and audits, conduction-driven viscosity and pressure, and late-observer refusal.


## One scenario control transaction for fluid inlets and outlets

`scenario_open_boundary::rebind_open_faces3` composes the existing velocity and
pressure preparations, then calls the native boundary owner once. Its six-slot
`ScenarioBoundary3` selection admits one control family per face; all selected
histories share one physical time and material scaling. Unselected faces keep
their complete existing controls. Selected spatial profiles retain each
original adapter's explicit refusal rule.

A late pressure failure, velocity-speed refusal, native kind/collision refusal
or cancellation cannot advance only one side of the flow. Both specialized
D3Q19 and generic D3Q19/D3Q27 owners retain their own collision, reconstruction,
force correction and step accounting. The transaction changes held targets;
it does not step the fluid, infer a geometry mapping, enforce total mass-flow
balance or claim a pressure-work certificate. Fixed six-slot preparation
copies no histories or population fields. G0/G3 regressions in
`tests/scenario_open_boundary.rs` exercise mixed native velocity/density
response and all-or-none refusal, including the original collision Mach bound.


## Hardness-dependent changing contact area in finite heat feedback

`friction_heat::plastic_wear::PlasticSlidingHeatContact::new_evolving_area`
opts the existing source-backed sliding/Archard/thermal owner into changes of the
positive-pressure cells on its original finite-gap grid. The complete projected
grid area must fit both uniform thermal face areas. Every nonlinear trial starts
from the exact accepted plastic, signed tangential and removed-volume histories;
recruitment and opening never reset or heal them. Original material cards still
supply temperature-dependent hardness, dry friction and caloric state separately.

Frictional dissipation and native Archard increments use the refreshed beginning
pressure. The final worn pressure footprint determines actual endpoint contact
area, source-derived half-cell conductance and ambient exposure through the same
backward-Euler contact-patch ingress as all other thermal links. Reports retain
the previous, refreshed and final areas and the exact accepted conductivity.
Successful candidates own their successor route and active mask, so the next
interval validates the accepted coverage and physical conductance. The default
constructor continues to require its original active cells.

This resolves total area between two uniform finite caloric receivers; it does
not resolve lateral thermal fields, constriction resistance, contact residence
or swept patches. The native beginning-pressure wear/endpoint-area split requires
physical time and grid refinement, and an activation discontinuity may refuse
nonlinear convergence. Complete detachment refuses an undefined massless
interface temperature. This base adapter still holds thermal mass and geometry
fixed within its explicit original-inventory removal envelope. Source regressions
exercise thermal-softening recruitment followed by wear-induced opening, actual
friction/volume/heat accounts, original irreversible histories, incomplete face
support and late source refusal using independent two-cell elastic/heat references.

## Runtime controls through pressure observation

`PressureRenderer` admits and applies batches of the existing `ControlDelta`
inputs at the current physical source boundary. Validation never advances time;
a refused batch leaves held inputs, pending events and control logs unchanged.
Producers without an explicit compatible input topology refuse nonempty batches.

Boxing and borrowing a producer preserve its input sink. `ScheduledRenderer`,
`EnsembleRenderer` and `DecimatedRenderer` forward controls to their original
retained mechanics. Decimation does not scale force, reindex source voices, move
source events, or reset the pressure filter. Both clocks remain unchanged by an
input update, and a poisoned producer refuses controls as well as rendering.

Live scheduled inputs apply immediately. A pending event on that same source
sample still runs in the schedule's existing stable order before the next sample,
so a scheduled assignment to the same input takes precedence. Live admission
reserves log capacity for the new batch and every remaining scheduled record;
the next pressure callback does not allocate merely because live logging consumed
the schedule's original capacity. This is a cold control operation, not a new
allocation-free control-thread guarantee. Schedule-only applied-event history
remains separate from the context's complete control log.

`ImpactPressureRenderer` accepts `SetModalForce` for its sole source slot zero.
The mode index, N/sqrt(kg) units and force ceiling are those of its original
`ImpactSource`, including a retained spatial thermal network. The complete batch
is checked before any held force changes; no mechanical, caloric, material,
radiation or propagation history is rebuilt. A force update uses the same held
input path as the concrete renderer's original `set_forces` operation.


## Exact reciprocal longitudinal piezoelectric storage

The longitudinal rod exposes the complete original quadratic storage block in
[extension, momentum, charge] order:
Q = [[k_D, 0, -theta/C], [0, 1/m, 0], [-theta/C, 0, 1/C]].
Its same source-bound energy and amplitude checks remain authoritative.
Native descriptor and prepared solvers can use analytic charge derivatives;
they do not need dimensionless finite-difference probes into microscopic
electrodes. The continuous storage tangent retains both reciprocal cross terms.
The quadratic declaration is exact for this fixed-state, small-signal rod and
is not a linearization of ferroelectric, hysteretic or nonlinear storage.

The G0 regression in tests/piezoelectric.rs compares every quadratic effort
with the original material response, independently reconstructs stored energy,
checks reciprocity and the analytic directional action, and verifies refusal
outside the original strain envelope. This adds no measured material
validation or hard real-time timing claim.


### Native circuit charge scales

The piezoelectric circuit now declares its complete reciprocal quadratic
differential block to native descriptor admission and stepping. Original
circuit storage, charge/strain cross terms, fixed-charge stiffness and inertia
are retained; algebraic multipliers carry zero storage. Analytic Jacobians
avoid generic finite differences of physical charge, including microscopic
electrodes. Material support and independent closed-work acceptance remain
authoritative. Regression definitions in tests/piezoelectric/circuit/charge_scale.rs
compare subfemtocoulomb response with independently eliminated loaded-rod
equations, both poling signs, zero coupling, domain refusal and retry.
These authored ceramic-scale fixtures are not measured material calibration.


### Prepared exact longitudinal current drive

`piezoelectric::prepared::PreparedPiezoelectricCurrent` prepares the original
three-coordinate [extension, momentum, free charge] system with the shared
`fs-phs::linear_exact::prepared::PreparedQuadraticFlow`. It retains every
reciprocal material entry and consistent rod inertia. Cold preparation has
explicit dimension, work and memory budgets; warm steps reuse fixed scratch
and the cached full flow and voltage integral. They perform no owner allocation,
Newton iteration, source clone or new matrix exponential.

Held force/current steps include the open-circuit case at zero current.
Mechanical supply is force times actual extension change; electrical supply
is current times the integrated trajectory voltage. The mean sensor voltage
is not an endpoint trapezoid. A computed conservative full-window energy bound
checks both original strain and field limits, including excursions between
equal endpoints. Independent charge/work acceptance and final cancellation
precede the sole physical state publication. The bound may refuse a physically
admissible step conservatively; it is not an interval-arithmetic certificate.

Synthetic G1/G4 definitions in tests/piezoelectric/prepared.rs compare an
independent held-current forced oscillator and its exact voltage integral,
repeated open-circuit frequency/charge, hidden interior excursions, and atomic
cancellation/retry. Exactness applies to the fixed linear held-input model
up to floating-point error; the existing spatial rod reduction, material
validity and absence of a hard real-time benchmark remain explicit limits.


### Moving-support spatial piezoelectric sensing and harvesting

`TetPiezoelectricTransient::step_with_support_motion_candidate` accepts one
endpoint velocity for every original constrained mechanical DOF. Support
coordinates retain their topology, while displacement follows the same midpoint
kinematics as the free coordinates. Full consistent-mass, constitutive and
reciprocal electrical coupling columns remain in the original solve. The
existing initializer still requires stationary supports initially.

Independent mechanical supply is applied-force work plus support-reaction
work. Returned midpoint support forces include inertia and intrinsic material
memory, with zero entries on free coordinates. Missing, duplicate, free,
nonfinite or unrepresentable support values refuse. Stationary APIs refuse
to silently stop a moving support; an explicit endpoint velocity accounts
for its deceleration and returned work.

The same motion is available through the resistive load owner's
`step_with_loaded_support_motion_candidate` and the reactive R/L/C owner's
`step_with_support_motion_candidate`. Original resistor heat, reactive storage,
field storage and source work participate in the enclosing balance.

The finite-temperature spatial owner exposes the same support motion through
`step_with_thermal_support_motion_candidate` and
`step_with_material_thermal_support_motion_candidate`. Every thermal proposal
starts from the accepted mechanical, electrical and material-memory history.
Actual support work enters joint supply once, while resistor and intrinsic
ceramic losses reach their explicitly bound finite caloric cells. The report
retains complete midpoint support reactions. Source-domain, inventory,
unchanged-coefficient, cancellation and physical-energy gates govern joint
publication. Without external resistors, nonlinear thermal trials reuse one
immutable native moving-support candidate.

The three synthetic journeys in
`tests/piezoelectric/spatial_transient_load/thermal/relaxation/support_motion.rs`
compare independent Maxwell/Debye, resistor-feedback, caloric and inertial
references through acceleration, continuation, reversal and stopping. Native
Rust compilation/execution remain unverified under the editing environment's
disk/memory limits; syntax/format and independent analytic checks were completed.
G1/G4 definitions in spatial_transient/support_motion.rs cover analytic rigid
inertia/reactions, independently condensed base-excited sensing, resistor and
capacitor harvesting, and transactional refusals. Small-strain assumptions,
source validity and temporal/spatial refinement obligations remain.

## Annular no-slip supply passages

`FixedChannel::Annular` carries the original `AnnularChannel` into the shared
two-sided and DG-backed aperture supply connection. Both physical radii, all
radial velocities, the inner and outer no-slip conductances, and the actual
gas source survive the connection. Pressure condensation, stored fluid energy,
physical viscous heat, numerical loss and gas-intervention work use the same
native momentum owner as the standalone annulus. No equivalent circular
radius or hydraulic-diameter resistance replaces the annular solution.

The common channel API also retains independent axial inner/outer wall drive;
coupled supply constructors retain their existing stationary-wall semantics.
Geometry is fixed and concentric, the fluid is incompressible and Newtonian,
and axial flow is fully developed. Rotating cylinders, normal wall motion,
eccentricity, entrance effects and transition are outside this model.

The existing DG supply startup/reversal/passive-release battery includes the
annular geometry. `tests/annular_channel.rs` exercises momentum-preserving
gas changes through the common channel API alongside the original analytic
Poiseuille, logarithmic Couette and narrow-gap checks.

## Second-order spatial moisture diffusion

`MoistureCartesian::advance_sdirk2` applies the existing Alexander two-stage
SDIRK tableau to the original three-dimensional finite-volume diffusion owner.
Both implicit stages resolve the pinned principal diffusivity laws and reuse
the native seven-point operator, checked CG, nonlinear controls and source
admission. Kirchhoff mode retains the original integrated D(M) face law,
including a degenerate dry state.

The first-stage contribution is evaluated from equal-and-opposite actual face
rates. The endpoint account contains (1-gamma)*dt times those rates plus
gamma*dt times the final rates through each of the six exterior faces.
The intermediate report is diagnostic and must not be added again. Dormant
surface transfer and impermeable directions cannot enlarge the admitted
second-stage moisture range. A nonphysical stage returns `StageRange` before
material evaluation; it never clips water or publishes a partial field.

G1 regressions specify anisotropic discrete-mode convergence, six distinct
boundary accounts, original nonlinear source refresh and analytic Kirchhoff
wetting. G0/G4 cases cover exact equilibrium, deterministic replay, dormant
reservoirs and cancellation. These are numerical checks with authored material
data, not experimental validation. Temperature and external faces are constant
within an interval. Ordinary `advance` remains backward Euler; no second-order
finite-atmosphere/solution coupling or global error bound is implied.

`MoistureCartesian::advance_sdirk2_adaptive_to` applies the same complete-duration
controller as the slab SDIRK path. The accepted two-half result uses the
second-order local-error factor (one third for equal halves), cubic step-growth
criterion and two-stage work charge per attempted endpoint. Range failures may
shorten the interval; source/domain failures remain terminal. Only accepted
RK-weighted six-face transfers contribute to the fixed whole-window water
budget. Rejected probes and accepted prefixes of a later failed call remain
private. Focused regressions specify the analytic error ratio, stiff-drying
recovery, complete-stage charging and late work-limit refusal.

## Simultaneous finite-wall, condensate and gas endpoints

`WallFilmThermalEnvironment::advance_implicit` couples original wall enthalpy,
volatile inventory and gas temperature in one shared Newton solve. The material
wall-film endpoint owns saturation, dry/wet complementarity, upstream vapor
enthalpy and finite-gas storage. The conduction endpoint owns wall calorics,
convection, radiation, source conductivity and the original Biot gate. Neither
an inner gas inversion nor an inner wall solve is performed during nonlinear
trials. After convergence each original caloric owner is reconstructed once,
and the complete phase/heat residual is checked again against both actual
returned states before the joint candidate is returned.

The finite wall receives condensation heat and supplies evaporation heat through
the same source-relative port as the material endpoint. Sensible gas-wall heat
has one reciprocal debit and credit; radiation remains an external enclosure.
Existing and newly formed film storage uses the liquid's original chart, never
the wall heat capacity. Exactly dry film has no thermal storage. At complete
evaporation, outgoing vapor still carries the wall-temperature-dependent
enthalpy. Wall melting/freezing plateaus remain primary-enthalpy states: latent
storage can change phase fraction at constant temperature without introducing
an apparent heat capacity or a duplicate heat source.

Numerical controls bound the shared Newton iterations and residual evaluations,
including the predictor, fresh converged trial and reconstruction rechecks.
Physical source, temperature, donor, mass and energy allowances remain those
of the original environment. Signed and absolute wall/film/gas balances retain
source-combination and boundary-port rounding; independent physical gates apply
after the numerical norm converges. Failure, source refusal, exhausted work or
cancellation returns no partial wall, film or gas state.

This opt-in endpoint path leaves the existing split advance and its continuation
semantics unchanged. It is first-order in physical time. The film exhaustion
time freezes the endpoint rate and requires timestep refinement for an accurate
continuous event time. Fully wetting, pure-liquid, fixed-area, fixed-gas-volume
and uniform-wall assumptions remain explicit. It does not add nucleation,
boiling, drainage, spatial gas gradients, displacement work or an entropy
theorem. Focused coupled regressions check actual finite-wall/gas feedback,
complete dry-out, retained latent wall storage and transactional refusal.

### Local discrete magnetic work in nonlinear voice coils

`speaker::magnetic` supplies a storage-owned symmetric magnetic chord for
`H_m = (lambda - psi_pm(x))^2 / (2 L(x))`. Endpoint-mean current and analytic
polynomial linkage/inductance secants close magnetic work in linkage and
excursion alone, including permanent-field and reluctance force. Independent
cone/panel inertia and reciprocal cavity energy keep their original quadratic
mean efforts and represented-midpoint energy remainder. The rule remains
well-defined at coincident displacements without dividing by their increment.

The original descriptor solver consumes the analytic endpoint derivative of
this complete rule. Constant `Bl`/`L` profiles retain the exact-quadratic
factorization. The existing thermal coupling observes the same discrete
current for coil-resistance heating. Source-support and arithmetic refusals
remain transactional; finite-step response still requires timestep refinement.
This displacement-dependent reversible law adds no current saturation,
hysteresis, skin/proximity loss or measured coefficient claims.

Focused regression definitions cover local and total work, independence from
unrelated mechanical momenta, reverse-step symmetry, endpoint derivatives,
coincident states, the constant-law path and invalid-support refusal.

## Annular DG terminal passages

Fixed DG terminal passages accept concentric axial annuli through
`DgTerminalChannelKind::Annular` and `dg_terminal_annular`. The original annular
momentum owner retains the actual local gas, unequal radial cell masses, both
no-slip walls, viscous history and separate pressure, kinetic, viscous, numerical
and material-change work. Admission requires positive ordered radii,
representable area, radial gap <= length/10, bounded radial cells and a resolved
speed Mach ceiling <= 0.3. Attached-area and combined-state budgets still apply.
Harmonic coupling retains every radial state. Finite-gas inventory includes the
actual annular passage volume, with main and branch gas ownership kept distinct.
This is a fixed concentric, axially developed Newtonian reduction. Eccentricity,
radial moving-pad openings, entrance jets and transition prediction remain
outside its scope. Existing terminal-group, gas-refresh, caloric-heat and
harmonic-field regressions now include annular geometry and independent ring
mass and logarithmic-shear calculations.

## Continuous simultaneous wall-film-gas windows

`WallFilmThermalEnvironment::advance_implicit_window` advances the coupled
enthalpy/volatile/gas endpoint over the complete requested duration. It shares
the original window clock, cumulative source apportionment and accepted-step
account with the split continuation. Wall internal heat and gas external heat
are full-window joules, assigned once by represented interval endpoints.
Accepted wall storage, liquid-film storage, gas internal energy, phase transfer,
sensible heat and radiation remain separate physical accounts. The implicit
gas storage already contains its sensible supply and is never counted twice.

`WallFilmImplicitWindowControls` bounds attempts and the complete window's
Newton iterations, residual evaluations and gas inversion work. Coarse or
rejected proposals cannot renew these budgets. Gas finalization reserves its
remaining original inversion bound; success refunds unused iterations, while
a failed finalizer retains the reservation when no partial receipt is available.
The report distinguishes known iterations from conservatively charged work.
Only an original phase donor `StepSizeLimit` may shorten and retry an interval;
source, caloric, nonlinear, balance and cancellation refusals remain terminal.

Every accepted step contributes to compensated whole-window balances. Final
admission independently checks actual wall storage, total and volatile mass,
carrier preservation, source allocation and accumulated absolute energy and
mass defects. Rejected proposals contribute no heat, mass, events or physical
time. A late failure publishes none of the accepted prefix. Exactly dry film
may form, evaporate completely and form again across successive windows;
original wall phase plateaus retain their enthalpy and phase-fraction storage.

This fixed-grid, first-order continuation inherits the endpoint model's
uniform-wall, pure-liquid, fixed-area and fixed-volume assumptions. Endpoint
rate exhaustion estimates still require physical timestep refinement.
Four focused groups cover manual endpoint composition and source allocation,
dry/wet/dry/wet continuation, native wall latent storage and late atomic
refusal; existing window and single-endpoint regression groups also pass.

### Second-order Cartesian diffusion with finite vapor

`FiniteMoistureAtmosphere::advance_cartesian_sdirk2` and
`advance_shared_cartesian_sdirk2` solve the actual finite vapor inventory at
both SDIRK stages. Shared specimens use one RH per stage, stable-ID ordering,
their original anisotropic/nonlinear material owners and one compensated sum
of local water changes. Stage two retains the ORIGINAL gas and solid storage
references and adds the explicitly weighted first-stage physical face fluxes;
stage-one uptake is never debited twice. Its first-stage transport receipt is
diagnostic. Endpoint face transfers, vapor, carrier species and gas energy
changes describe only the original-to-final physical interval.

Sorption and optional surface emission refresh at each trial RH. The accepted
stage envelope includes original field extrema and only equilibria connected
by nonzero first-stage or current second-stage physical face coefficients.
Both pre-source and final admission use that envelope. Inadmissible stages
refuse without clipping or committing a partial atmosphere/specimen state.
Every scalar-bracket trial must remain inside admitted source support; an
inadmissible extreme trial may require a shorter timestep.

The original prescribed-temperature, fixed-volume dilute-vapor scope remains.
There is no condensation, sorption heat, moving interface, global temporal
error certificate or experimental-material validation implied here. G1
manufactured coupled relaxation, nonlinear flux/source and shared-inventory
fixtures check the method; G0/G4 fixtures cover equilibrium, cancellation and
stiff-stage refusal. Existing backward-Euler entry points retain their method.

The finite-vapor Cartesian SDIRK path also provides
`advance_cartesian_sdirk2_adaptive_to` for a complete requested duration.
The existing controller checks moisture, RH and vapor mass independently with
second-order local step-doubling scaling, charges both implicit stages per
proposal, and sums only accepted physical half-interval RK transfers. Final
water, carrier and total-mass checks admit both inventories together. Failed
stages, rejected probes and a late work-limit failure publish no partial
duration. G1/G4 fixtures independently check the coupled relaxation, all three
error coordinates, stiff-stage retries and complete-stage work reservation.

## Adaptive time accuracy through wall and film latent storage

`WallFilmThermalEnvironment::advance_implicit_adaptive` compares one simultaneous
endpoint with two consecutive half endpoints from the same accepted state.
The original thermal first-order controller supplies the local-error correction,
including unequal represented halves, and bounded step growth/reduction.
Only the fine pair enters the shared physical window account and retained
history; all three solves and rejected comparisons spend the same global
Newton, residual, gas-inversion and attempt budgets.

The componentwise error estimate observes original wall enthalpy, liquid phase
fraction, wall and gas temperatures, actual film inventory, every catalog gas
species, and liquid-film specific enthalpy when both endpoints are wet.
Relative allowances scale changes from the interval's accepted start, so an
arbitrary caloric datum cannot loosen accuracy. A latent wall plateau can
therefore require refinement even when its temperature is exactly constant.
A dry film has no invented caloric coordinate or positive mass seed.

Uniform body and gas source totals are apportioned at absolute cumulative
window endpoints. Coarse and rejected sources never reach accepted balances.
Only an original phase donor step-limit refusal permits admission subdivision;
other source, nonlinear, conservation and cancellation errors remain terminal.
Represented admission spans never exceed the source's maximum. A bounded
comparison/reduction count and minimum span prevent endless accuracy retries.

The report exposes accepted/discarded work, actual fine samples, component
differences, maximum accepted normalized error and an explicit next-span
suggestion. This remains a first-order physical discretization with a local
error estimate, not a global trajectory or continuous dry-out-time bound.
Known changes of forcing require a new authored window. The endpoint's
uniform-wall, pure-liquid, fixed-area and fixed-volume limits remain.

Focused regressions compare tightened controls against independent continuous
wall/film/gas balance equations, isolate phase-fraction refinement at constant
wall temperature, and verify fine-only source composition, complete discarded
work charging, bounded rejection and atomic refusal.

### Reciprocal acoustic guide loading of spatial piezoelectric transducers

`harmonic::acoustic::PlaneWaveAcousticPort` prepares an original exterior-face
pressure/flow map. The exact P1 outward normal weights b give U=b^T v and
traction -b p. One original `FluidAcoustics` owner supplies the coherent
density, adiabatic sound speed and bulk modulus at the evaluated state.
The retained uniform guide mode has Z=rho*c/A and p=Z U+2 p_in under the
original exp(-i omega t) convention. The guide section can include a declared
rigid area; moving-face geometry remains the original tetrahedral trace.

`solve_harmonic_with_acoustic_ports` adds -i omega Z b b^T and -2 b p_in to
the same sparse electroelastic pencil, including intrinsic relaxation and
external electrical shunts. Upper-triangle entries are mirrored exactly.
Pressure recovery retains complete fluid loading in support reactions.
Independent incoming/outgoing wave powers and signed p/U work remain distinct
from conductor and ceramic heat. Net outward acoustic power enters the same
physical work acceptance gate; incident sound may supply a real sensor/shunt.

Ports retain their exact model, unique consistently outward planar exterior
faces, complete guide area and explicit surface-coordinate work bound.
Duplicate/foreign ports, inconsistent geometry, invalid original fluid states,
out-of-band wavenumbers, excessive pressure/flow amplitudes and cancellation
refuse before publishing a response. Material evaluation occurs once per port
before sparse assembly, never in Krylov iterations.

This is a declared single uniform-mode, lossless semi-infinite guide model.
The caller owns its transverse-mode cutoff/accuracy domain and rigid remainder.
It is not free-space diffraction, a higher-mode/evanescent approximation,
viscothermal guide loss, a nonlinear acoustic model or empirical validation.
Four synthetic G1/G3/G4 regression definitions specify electrical actuation,
incoming-sound sensing with shunt heat, rotated pressure traces, rigid reflection
and refusal. Rust syntax/format and source/equation review passed; independent
scalar work calculations agreed. Native compile/tests remain unverified under
the editing environment's disk and memory limits.

## Source-bound saturating core, gap and leakage inductors

`material_circuit::inductance::MaterialReversibleInductor` maps an admitted
reversible H(B) source into the native circuit's bounded I(lambda) storage.
The declared uniform-flux reduction is Phi=A*B,
N*I=l*H(B)+R_gap*Phi, lambda=N*Phi+L_leak*I. Original source breakpoints are
transformed exactly, with both polarities of the isotropic zero-remanence law.
No sampled fit, constant apparent inductance or second time integrator replaces
that relation. Source snapshot, query context and fixed geometry stay inspectable.

Its stored energy is the same potential in both descriptions:
A*l*integral H dB + R_gap*Phi^2/2 + L_leak*I^2/2.
Operating points expose core B/H, current, differential inductance and separate
core/gap/leakage energies. These partition the native storage; they are never
additional energy debits. B is recovered from the retained transformed segments,
avoiding cancellation of lambda-L_leak*I when leakage dominates linkage.
The native right-linkage derivative owns source kinks, except at the final
endpoint where the left derivative applies.

`ReversibleInductorCircuit` binds these sources to original inductor branches,
sets their zero-field differential coefficients, and delegates complete
Kirchhoff assembly, initial consistency, analytic actions and audited stepping
to the existing `CircuitDae`. Nonlinear capacitors, reciprocal capacitance
blocks, other authored scalar inductors and disjoint linear mutual links keep
their native owners. No extra accepted history or time integrator is introduced.

The source must be a field-only reversible potential. Explicit thermal scaling
or additional response axes refuse, because freezing a thermal free-energy
snapshot cannot supply its magnetic entropy account. Source-domain exits,
duplicate ownership, scalar nonlinear/mutual overlap and unrepresentable
geometry refuse. Gap reluctance and leakage are explicit fixed inputs; this
reduction does not infer fringing, distributed field resolution, hysteresis,
eddy loss or temperature evolution. Graph resistance retains physical heat.

Focused regressions cover every source knot and polarity, independent
core/gap/leakage work, turns scaling, leakage-dominated observations, a native
saturating RL trajectory and physical/source refusals. Independent arithmetic
checks over 10,000 varied geometries agree with the source potential within
1.5e-15 relative error. Rustfmt parsed the implementation and tests; complete
source-to-native circuit tests remain unexecuted in the dependency environment.


### Prepared source-bound piezoelectric circuits

`material_circuit::piezoelectric::prepared::PreparedPiezoelectricCircuit` caches the original fixed-material index-one Kirchhoff system for a fixed positive clock. Native algebraic elimination retains every differential flux, charge and mechanical coordinate. Native exact quadratic propagation and independently prepared physical source observations preserve resistor feedthrough, separate electrical/mechanical supply and per-resistor heat. Original endpoint Kirchhoff rows, source-domain limits, material conductor power and total stored-energy balance are checked before transactional publication. Warm stepping performs no owner allocation, factorization, exponential or Newton solve. Singular algebraic blocks refuse explicitly without reducing conserved state; ordinary descriptor stepping remains available. This is the linear held-input model up to floating-point error, with a conservative numerical whole-window domain bound, not measured material validation or a hard real-time guarantee.


### Prepared fixed-voltage piezoelectric motion

`piezoelectric::prepared::PreparedPiezoelectricVoltage` owns the original rod and complete accepted extension, momentum and free charge under the terminal voltage inherited at preparation. Native exact quadratic propagation eliminates only the algebraic voltage constraint, retaining the short-circuit stiffness, consistent inertia and reciprocal force theta*V. Each candidate reconstructs charge from the original preparation anchor and independently checks the original full stored energy, force work, voltage-source work, induced charge and constitutive voltage before publication. Native whole-window mechanical norm admission bounds strain, and the fixed voltage with its explicit tolerance bounds electric field. Warm stepping uses fixed arrays and performs no owner allocation, factorization, exponential or Newton solve. Refusal and cancellation preserve the complete accepted body. This is the held-input linear model up to floating-point error with conservative numerical domain admission; changing to an incompatible ideal voltage requires another physical boundary model and cannot silently reset charge.

## Regional coexistence in the conserved finite water cell

`EquilibriumWaterCell::from_isochoric_state` admits an original prepared
`IsochoricWaterState` directly. Its mass, fixed specific volume, raw internal
energy, source mode and dense-root controls remain the source of the initial
M,U,V inventory. Initialization requires the exact bounded source controls
after applying the cell's absolute budgets per unit mass.

The original stream, finite-wall heat and sealed moving-volume paths now
preserve that source mode when rebinding V/M. Regional coexistence can therefore
evolve liquid/vapor masses, pressure and density above 623.15 K through the
same conservative ledgers and reciprocal pressure-work solve. Phase volumes
use the retained original regional endpoints; mass quality is never substituted
for vapor volume fraction. The constructor does not create a second EOS or
replace conserved energy with an approximate recovered source image.

Both native water-cell acoustic observers support the regional pair and use
exactly the retained thermal source density, with the selected original bulk
modulus. They expose the signed arithmetic density difference separately.
The equilibrium observer preserves the original Gibbs accessor and adds a
regional tangent accessor retaining its actual root/pressure defects. These
are explicit asymptotic acoustic closures; the observer does not add phase
kinetics or advance the thermal inventory.

Regional operation remains a bounded coexistence window below the critical
temperature, on one side of the 623.15 K source-fit seam. Losing either phase,
crossing that seam, or exceeding the original source/residual budgets refuses
the proposed interval. No noncondensable gas, slip, nucleation, moving interface
or thermal-gradient model is inferred. Stream transfer keeps its existing
held-donor approximation; moving-volume work keeps its existing trapezoidal
time accuracy. Coupled material and source errors remain distinct from
conserved mass/energy residuals.

Five added native integration tests cover 640 K heat/phase/acoustic updates,
two-cell source-retaining transfer, reciprocal isothermal phase expansion and
return using raw Delta U plus p Delta V, source-window/late-cancellation
rollback, and finite conducting-wall exchange on the same regional inventory.

All 23 exact native water-cell tests passed, including the five additions and
18 existing regressions. The source closure uses original feedback helpers and
genuine native material, duct and thermal dependencies; no whole-workspace
Cargo/DSR pass is implied.

## Finite-rate water compression in reciprocal cavity networks

`cavity_phs::relaxation::RelaxingWaterCavity` retains the original
`RelaxingRegionalWater` and `HelmholtzCavity` geometry. Native neck inertia uses
the exact thermal density and declared end corrections. Original equilibrium
bulk storage and one positive compression arm enter the existing
`RelaxationBranch` and `PreparedExactZoh` owners. Pressure and volume flow remain
reciprocal, including through native mechanical transformers.

Fixed-volume net inflow gives first-order excess mass and density. Pressure,
neck flow, storage and instantaneous compression/neck losses come from the
accepted native state. Relaxed or explicit frozen memory is supported. Native
prepared stepping retains integrated input work and dissipation with reusable
warm scratch.

This homogeneous small-signal lumped model retains caller-established
relaxation kinetics, explicit neck resistance and end geometry. It does not
evolve nonlinear phase inventory or temperature or infer exterior radiation
or viscosity. Four focused native regressions pass: original 100 kPa/20 MPa
density and inertia; independent refined mass/momentum/memory ODE plus
work/heat; actual mechanical-transformer power and total energy; native
admission and tiny density perturbations. Exact-source tests used genuine
material and pHS owners; no full workspace Cargo/DSR pass is claimed.

## Motor torque in source-derived directional sphere contact

`friction_heat::directional::sphere_plane::advance_source_sphere_plane_with_torque`
threads a prescribed constant world couple through the original native sphere
trajectory and its original source-resolved support and impact laws. Actual
torque-driven spin reaches preimpact speed queries; supported motor rolling,
breakaway and sliding reach the motion-consistent speed/path source closure.
The full static and kinetic material claims, original selection receipts and
plane-owned directional frame remain authoritative. The existing entry
delegates zero couple.

The returned native motion receipt separates signed motor work and angular
impulse from centre-force work, moving-plane work, sustained slip heat and
unallocated impact loss. Query callbacks see the latest accepted private
mechanics/history. A late source refusal or cancellation discards the entire
candidate, including earlier motor motion and collision source updates.

Four G1/G4 regressions in
`tests/rigid_directional_heat/sphere_plane/impacts/torque.rs` cover sourced
static rolling with zero path/heat, motor breakaway with actual accumulated
slip, torque-dependent preimpact source coefficients and bounce work, and late
source/cancellation rollback. The fixture ledgers were independently derived
and all changed Rust sources parse with native rustfmt. These tests have not
been executed against the current fs-mbd/fs-couple dependency closure; no
Cargo/DSR gate is claimed. This remains an isotropic rigid point-contact model
with a held world couple, explicit temperature/pressure owners, and no inferred
contact area or impact-to-heat allocation.

## Adaptive second-order diffusion with shared finite vapor

`FiniteMoistureAtmosphere::advance_shared_cartesian_sdirk2_adaptive_to`
extends the same second-order controller to separate native Cartesian specimens
sharing one vapor inventory. Each field retains independent local error and
six-face water accounts, including simultaneous uptake/desorption with zero net
reservoir change; RH and vapor mass are checked separately. Only accepted
RK-weighted half intervals contribute transfers. Both implicit stages count per
coupled probe, stable IDs determine output order, and final whole-window water,
carrier and total-mass checks admit all endpoints together. The local error
estimate remains asymptotic, without a global-error guarantee.

G1/G4 regression definitions cover coupled mean/difference modes, stationary
vapor with opposing transfers, ID permutations and late stage-budget failure.

## Coupled orifice law on the first N open tone holes (`open_hole_orifice N`)

The optional reed-card record `open_hole_orifice N` (N = 1..=3; after `constriction_noise`, before `compile_limits`; once
only) solves the quasi-steady sharp-edged orifice law `delta_p_i = K_i |U_i| U_i` (`K = rho / (2 C_vc^2 A_h^2)`, the
first-hole shunt's own law and coefficient) of the first N open tone holes TOGETHER instead of the first alone
(`RegionalDuct::with_open_hole_orifice_ports`; `RegisterHoleOrificePorts` in `regional/outlet/orifice_ports.rs`). The
holes couple through the bore, so the per-sample solve is an N-vector one, a damped Newton on `U = u_th - A Phi(U)` with
`A = R^-1 (I - D)(I + D)^-1`, `D` the direct taps of the N x N reflection kernels `S = (I + R Y)^-1 (I - R Y)` and `R`
the holes' numerical port references (the ordering `R Y`, not `Y R`, matters as soon as the references differ). `Y` is
the admittance of the linear TMM chain seen at the N orifice terminals
(`fs_duct::exterior::input_impedance_holes_and_ports`: a nodal ladder of the blocks' own admittances, a pure-shunt
block such as a sealed hole without series length sitting on its node; every other hole stays linear), and it equals the
inverse of the one-port network impedance for N = 1. The reed receives the sum of the N `s12` kernels applied to the
accepted drops; radiated hole sound still comes from the linear network, as for the first-hole shunt.

An absent row or `N = 1` runs the unchanged first-hole shunt, bit for bit (`RegisterHoleShunt::Single`), and an N above
the number of open holes uses all of them. For N >= 2 the realized matrix is made passive by bounding the matrix
itself: the direct matrix is clamped to power-wave norm 0.999 and the whole tail scaled by the largest factor (found by
bisection, 2^-14 resolution) for which `max_w || D + lambda T(w) ||` stays at most 0.999 on the 2n grid. The
first-hole line's own rule, `|D| + max_w |T(w)| <= 0.999`, ignores the phase between `D` and `T`: it costs 0.1 % of the
tail for one port but 0.4 % for two and 3.3 % for three on the test bore; with it the coupled flows departed from the TMM
chain by up to 15 % at 100 Hz on a 4096 grid, with the bisected rule by at most 3.4 %. That legacy rule is kept (`PassivityRule::Triangle`) only so the N = 1 equivalence
test is exact; N = 1 in the card never reaches it. The row needs one constant-medium run: a `gas_field` or molecular
bore refuses by name at parse (N > 1) and a multi-region duct refuses at build. Cost: `N^2 + 2N` kernels of the grid
length instead of 3, and one more frequency sweep at build.

Evidence (G0/G1): `fs-duct` `compose::port_tests` (the ladder's admittance is the inverse of the one-port network
impedance to 1e-9, reciprocal, passive, independent of the port count, and the sign of its off-diagonal term fixed by a
coincident-pair test); `orifice_ports::tests` (the solver's residual over 600 random passive systems and a stiff drive;
the realized matrix is passive and is the admittance scattering, not its transpose; the time-domain steady state with
linear drops equals the network built from the realized kernels to 1e-14 and the TMM chain to a few percent at a 4096
grid; the general solver run as one port follows the first-hole line through a strongly nonlinear drive to 1e-6);
card level `reed::orifice_tests` (an absent row and `open_hole_orifice 1` are bit-identical in audio and reed state, a
one-hole duct ignores the count, 2 and 3 differ from 1 and from each other, bad counts refuse).

## Second-order diffusion coupled to finite solution activity

`FiniteMoistureSolution::advance_cartesian_sdirk2` and
`advance_shared_cartesian_sdirk2` solve native Cartesian moisture fields and
the same finite solvent inventory at both implicit stages. Actual water uptake
concentrates the solution; desorption dilutes it. Every trial refreshes solution
activity, original sorption and surface transfer, and nonlinear principal
diffusivity. Original-to-final compensated water increments close the reservoir
once. Stage-one fluxes enter the second-stage load exactly once; the retained
first-stage receipts are diagnostic, never additional solvent debits.

Dissolved ideal/nonideal solutions, neutral equilibrium solubility and ionic
equilibrium keep their original thermodynamic and total-component owners.
Second-order accuracy assumes a smooth constitutive/phase branch. Kinetic
solubility is explicitly refused because its endpoint-frozen phase update does
not supply the required coupled RK equations. Equal prescribed temperature and
pressure, negligible vapor storage and the original source domain remain in
force; no liquid-immersion permeability, ion transport, thermal evolution or
mechanical feedback is inferred. Shared participants retain unique stable IDs,
native geometry, source laws and aggregate cell admission. Borrowed original
fields and loads avoid cloning either for each scalar reservoir trial.

`advance_cartesian_sdirk2_adaptive_to` completes one requested duration with
the original temporal controller. Native content, solvent mass and activity
receive independent error scales; both stages spend work even in rejected
probes. Only accepted physical half-interval transfers enter the face ledger.
One fixed whole-duration water budget admits both endpoints together; a failed
duration returns no partial bath or field.

G1/G4 regression definitions cover independently integrated finite-bath
relaxation/refinement, equilibrium phase buffering and solute conservation,
opposing shared specimen transfers, kinetic/cancellation refusals, and adaptive
whole-duration accuracy and work limits.

## Condensate films on the original spatial thermal graph

`finite_gas_thermal::wall_film::deposition::advance_deposition_with_wall_film`
attaches one independent Dry/Wet film to an existing native thermal cell.
The original `BackwardEulerStep` supplies every cell's caloric storage, contact
flux, source and prescribed bath. The original material wall-film endpoint
supplies evaporation/condensation, zero-inventory complementarity, liquid storage,
vapor-carried energy and finite-gas feedback. A shared Newton solve uses native
cell-energy increments, one volatile inventory coordinate and gas temperature.

The selected patch occupies actual exposed wall area. Its gas convection is
reciprocal; the remaining area retains the original ambient cooling. Native
reservoir coverage and contacts remain authoritative. An ambiguous callback
replacement of that occupied ambient channel refuses. Other endpoint sources,
link conductances, internal transfers and bath replacements remain available,
and their callback sees the actual trial wall graph and gas.

The underlying wall keeps its mass, material identity and original phase chart.
Liquid-film h(T) comes from its separate source; melting/freezing storage in the
wall is never treated as film mass or apparent heat capacity. A conductive
neighbor may supply evaporation heat or absorb condensation heat while the
wall stays on its latent plateau.

No gas or wall caloric root is nested inside residual evaluations. After one
original gas reconstruction, a fresh native heat/phase evaluation checks the
actual joint endpoint. The independent external-heat audit runs only on this
final replay, reusing the same resolved coefficients without another callback.
Local cell defects, the complete signed/absolute energy account, source-port
roundoff and original species/mass balances all gate publication. The shared
implicit work meter also supports enclosing continuations.

This event-free first-order endpoint preserves native exact-target event
semantics and fixed geometry between events. Each cell remains a uniform
control volume requiring spatial/time refinement; the film is fully wetting,
pure liquid and thermally equilibrated with its wall patch. It adds no drainage,
moving contact line, boiling or liquid-mixture activity model. Focused tests
cover analytic conductive supply on a wall fusion plateau, dry condensation,
the uniform-wall limit, manufactured partial-area/source feedback and atomic
physical/source/event/work refusal.

## Native thermal events with a continuing condensate film

`deposition::timeline::advance_deposition_with_wall_film_timeline` advances
the spatial wall-film endpoint through the original thermal schedule. Births,
contacts, separations and boundary changes apply in authored equal-time order,
including events at the initial and requested final time. Every heat/phase
interval stops before a native event; the event pass permits zero thermal
substeps, so it cannot insert an independent sensible-only heat advance.

The initially active wall keeps the same physical wetting patch and separate
film inventory. New contacts must leave enough native exposed area for that
patch; incompatible coverage refuses the whole candidate instead of hiding,
draining or discarding condensate. Unselected cells may be born during the
window and subsequently exchange heat through original contacts. Birth mass
and energy retain each new cell's original material datum. Final wall storage
subtracts each cell's initial or birth enthalpy before multiplication and
summation, so incoming formation energy cannot masquerade as latent heat.

External gas heat is a full-window uniform source, apportioned by cumulative
physical endpoints. Wall sources retain their native watts and boundary events.
Only original phase step-limit refusals shorten a trial; other source, nonlinear
and conservation failures propagate. All proposals share the original implicit
work meter, the complete-call attempt limit and the native thermal substep
allowance. Accepted intervals alone contribute physical accounts. A final
source, coverage, conservation or cancellation failure returns no partial graph,
gas, film or event history.

Compensated totals retain accepted phase, sensible, storage and external-source
accounts. Final admission checks actual native wall storage against interval
sums, exact total gas-source allocation, film/volatile/whole-gas mass, carrier
preservation, and cumulative absolute mass/energy defects. The final payload
belongs to the last accepted heat endpoint before any final-time events.

This bounded continuation remains first-order with a fixed maximum interval;
native event clipping and donor admission do not estimate temporal accuracy.
The three focused groups compare manual event-cut composition, isolate a
large independent birth datum from actual conductive/latent heat, and exercise
late physical/source/work/cancellation refusals.

## Adaptive second-order diffusion with shared finite solutions

`FiniteMoistureSolution::advance_shared_cartesian_sdirk2_adaptive_to` completes
a requested duration through the existing shared temporal controller. Every
original specimen, solution activity and solvent inventory receives separate
local error control; opposing transfers cannot cancel specimen error estimates.
Complete-window accounts include only accepted RK-weighted transfers, preserve
the native nonvolatile stocks, and charge both implicit stages of every proposal.
All fields and the one finite reservoir are admitted together after water
conservation checks. Native kinetic-phase refusal remains; temporal estimates
assume smooth source and algebraic phase branches and certify neither global
error nor phase-crossing order.

G1/G4 regression definitions cover opposing fields with a stationary pure bath,
finite concentration/activity feedback against an independent continuous
reference, accepted signed transfers and late stage-work exhaustion.

## Adaptive exponential binary diffusion in the live material assembly

`MaterialExchangeState::adaptive_binary_diffusion_step_candidate` and
`advance_with_adaptive_binary_diffusion` connect the original finite binary
gas network to the existing reaction, thermal, electrical and piston owners.
The gas stage refreshes actual source states and partial enthalpies, resolves
local species/temperature error, and enforces separate trial/pair-work limits.
Every accepted gross transfer marks both chambers' closed two-port audits open,
including counterdiffusion with zero net mass.

The original assembly balance encloses the pre-diffusion inventories. Actual
diffused species feed the original reaction owner once; transported and reaction
enthalpy remain in their gas inventories and are not redeposited as external
heat. The common feedback and whole-window mass/energy gates must pass before
gas, thermal field, charge, motion and clock publish together. Source failures,
pair-work exhaustion, late feedback refusal and cancellation preserve the
entire accepted assembly.

The complete diffusion/reaction/heat-motion split remains first order, and the
gas-only local estimate does not certify coupled splitting error. This path
retains the original binary molar-frame and fixed-volume transport substage;
pressure/thermal diffusion, resolved nozzle flow and a simultaneous
moving-volume constitutive solve require their existing distinct owners.
G1/G4 regression definitions cover an independent stiff exponential limit,
original-owner composition with the full feedback solve, reaction consumption,
zero-net-mass open audits and transactional late failures.


### Continuous prepared piezoelectric force/current tones

`PreparedPiezoelectricTones` advances continuous fixed-frequency force/current programs through the original three-state source-bound rod and native extended quadratic-flow owner. One shared numerical cosine/sine carrier per frequency retains both ports' independent phase offsets without introducing physical storage. Native signed quadratic observations integrate mechanical and electrical work independently; native linear observations integrate current and terminal voltage. The actual accepted carrier norms enter a conservative whole-window strain/field bound, including numerical phase-norm drift. Body, phase, accepted clock and reports publish atomically after original source-domain, charge, energy and cancellation checks. Warm stepping allocates no owner buffers and performs no exponential, factorization or Newton solve. Scope is the fixed-material linear sinusoidal/DC source model up to floating-point error, without phase normalization, energy rescaling, empirical material validation or a real-time deadline claim. Six focused regression definitions cover resonance, weak coupling, multitone interference, independent source work, DC parity, cancellation and whole-window/budget refusals; native execution remains a separate verification obligation.


### Transient piezoelectric acoustic guide coupling

The original spatial midpoint solve accepts real incident/outgoing guide waves through the same `PlaneWaveAcousticPort` pressure/flow trace used by harmonic acoustics. The actual immutable `FluidAcoustics` owner supplies density and adiabatic sound speed. Reciprocal acoustic traction, consistent inertia, charge constraints, optional Maxwell/Debye ceramic memories and physical resistive harvesting share one native solve; prescribed support motion pays its complete reaction work. Acoustic outflow is audited against independently recovered incoming/outgoing wave energies and is counted once, separately from ceramic and resistor heat. Time-domain midpoint energies use no peak-phasor factor of one half. The accepted endpoint retains constitutive forces without spurious acoustic history. Foreign/overlapping traces, unrepresentable waves, amplitude/domain violations and cancellation refuse transactionally. This is a fixed-geometry uniform-guide model with caller-declared waveform/interface bandwidth; a timestep is not a spectral certificate. It does not imply free-field radiation, diffraction, a spatial acoustic volume, or time-discretization convergence. Three focused synthetic regression definitions cover sound-to-electric conversion with/without intrinsic memories, moving/rotated interfaces, reflection and domain/cancellation behavior; native execution remains a separate verification obligation.


### Finite-temperature full-tensor spatial piezoelectric dynamics

`spatial::transient::load::thermal::tensor::ThermalTensorPiezoelectricBody` couples the full native 3-D spatial fields and consistent-mass motion to one uniform finite ceramic temperature, with original per-cell tensor orientation and complete source-bound caloric inventory. Temperature proposals refresh the complete affine stress-charge C/e/epsilon source tensors and reassemble the same native operators under original budgets. The model may own refreshed tensors while borrowing unchanged geometry; identical-temperature trials reuse accepted operators. Each changed-temperature trial first solves electrostatics at fixed original displacement/free charge, restores physical momentum/history, and then runs the original shunt/current/support-motion midpoint solve. Signed reversible field/caloric transfer is derived from the full Helmholtz and entropy endpoints; actual resistor Joule heat is deposited only through its original receiver. Independent acceptance checks total field internal energy, kinetic energy, all finite caloric energy and actual external work/heat. Complete source/context/density/inventory checks, declared strain/electric-displacement bounds and a whole-interval caloric-plus-field capacity floor precede coupled Newton publication. This is fixed-reference geometry with one uniform ceramic temperature, fixed zero-voltage grounds and charge/current electrodes; intrinsic memory, nonzero ideal-voltage constraints, expansion, pyroelectric offsets, phase transitions and independent cell temperatures require additional accounted models. Endpoint-temperature coefficients make the thermal feedback generally first order in time. Four focused synthetic regression definitions cover moving/rotated supports, reversible caloric changes, heating-dependent capacitance at conserved charge across mesh refinements, separate finite shunt heating and transactional refusals; native execution remains a separate verification obligation.
