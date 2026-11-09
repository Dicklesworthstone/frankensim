<!-- fs-couple contract, part 8 of 8; a higher number is newer. Index: ../CONTRACT.md. -->
## Driven full-curve Hertz friction with original thermal sources

`friction_heat::hertz::rate_dependent::driven` connects the original local-speed
Hertz initial-value trajectory to all selected `SourceThermalFriction` curve
parameters. Each trial refreshes the same pinned claims at the solved common-face
temperature and advances the captured body pair, held external wrenches and
explicit contact mode. The existing finite-cell phase solve owns physical heat
and source admission; motor work remains a separate mechanical supply.

The original mode-aware full/two-half marcher carries accepted bodies, mode,
caloric state and source receipts together. Native estimated travel and heat
errors enter the existing local comparisons alongside whole-span independent
mechanical, thermal and combined work balances. No parallel source model, thermal
store or adaptive controller is introduced. Coarse/fine mode disagreement,
source-domain failure, unresolved original-law events and late cancellation
refuse without publishing any private candidate.

Spatial pressure, mechanical integration and thermal-window errors remain
independent axes. Geometry, normal load, radius, inertia, footprint and external
drive stay held per phase. Genuine Coulomb retains its original event path and
explicitly unavailable adaptive error estimates; other curves refuse unresolved
driven zero-rate events. No automatic weakening-root selection or general
nonconstant driven arrest is claimed. Direct tests compare exact affine motor
and finite-cell roots, an independent Stribeck spin pressure integral, continuous
feedback in both dry strength and speed slope with separate power integrals,
original source identities, static cooling and atomic retry.

## Driven Hertz friction, finite-cell heat and retained contact mode

`friction_heat::hertz::driven` proposes original two-body driven shear/spin
phases and adaptive finite-cell thermal continuation. Each common-face
temperature trial resolves the same original static/kinetic source claims
and restarts the captured body pair, external world forces/body torques, and
explicit resting/sliding mode. Only original physical friction loss heats the
finite cells. Signed full-body external work remains a separate supply,
including common-body acceleration while a contact sticks without heating.

The existing full/two-half thermal march carries the fine accepted mode,
original body momenta and caloric state together. Coarse/fine disagreement in
retained mode forces refinement even when numerical rate differences are small.
Drive work has a separate coarse/fine defect and compensated accumulated account.
Whole-span mechanical, thermal and combined balances separately check
`Delta K + Q - W_drive`, `Delta U - Q - Q_external`, and
`Delta K + Delta U - W_drive - Q_external`. Late refusal or cancellation
publishes none of the private body, mode, material or thermal candidates.

Geometry, inertia, load, radius and footprint remain held; native admission
bounds closed normal motion over the complete phase. Thermal source sampling
is per phase, with independently refined first-order curved mechanics. Neither
exact thermal threshold timing, driven rate-weakening branch selection nor
normal evolution is claimed. Focused tests cover static common-motion work,
powered spin/arrest, thermal breakaway and discrete-mode disagreement,
independent coupled motion/power/caloric refinement, and atomic retry.

## Scenario temperature instruments on accepted finite-cell heat

`scenario_sensor::TemperatureSensorObserver` connects original compiled scenario
instruments to `DepositionSimulation` through its existing `DepositionObserver`
contract. An explicit entity-labelled state-slot map selects distinct active
persistent cell IDs. Quantity, entity and dimension admission prevent pressure
rows, unrelated entities, duplicate cells and map-order changes from silently
becoming temperature probes. The caller still owns geometric correspondence and
patch weights; this adapter does not infer a point location or camera footprint.

Every native accepted substep supplies its original caloric-law temperatures,
physical duration and endpoint clock. The scenario's shared mount operator and
declared instantaneous/first-order response evolve on an affine reconstruction
of those temperatures. The same original law supplies temperatures through
latent-heat plateaus. Original boundary/contact events preserve sensor memory.
The observer adds no heat, conductance or time restriction and composes with
other physical owners through `CoupledObservers` on the same thermal timeline.

Initial response is explicit. Stale clocks or starting temperatures, invalid
support, material refusal and cancellation abort the existing thermal candidate
and its cloned instrument history. Cell sampling and field evaluation poll at
most every 256 entries; event checks retain the original bounded synchronous
event-hook semantics. Storage is linear in the admitted sensor state dimension,
with no growing history log. This is a passive observation of finite-cell heat,
not thermal sensor back-reaction, continuum sampling accuracy, stochastic noise,
calibration validation or a transient assimilation of uncertain sensor memory.
`observe_endpoint_candidate` also consumes an original implicit or multiphysics
thermal candidate without rerunning heat. It shares the same response owner and
requires the matching accepted start state and original duration. The caller
accepts both candidates together and explicitly assumes an affine temperature
path across that one interval; hidden substeps or events are not reconstructed.
`tests/scenario_sensor.rs` exercises real heating, mount bias, lag, patch probes,
source-event timing, native-clock rounding, restart, phase plateaus and rollback.

## Canonical magnetic scenario boundaries

`scenario_magnetic::magnetic_flux_boundary_at` and `solve_magnetic_flux3`
bind typed normal-flux histories to the original `SpatialMagnetostatics`
exterior faces. Inputs are uniform B.n densities in tesla with odd reflection
parity and an explicit inward/outward sign convention. Sampling preserves
physical-time frame validation, source policy and canonical quantity metadata.
The native owner retains actual face integration, per-component flux balance,
deterministic gauges, material sources, residual gates and stored energy.

Unsupported metadata, overlapping/interior faces, invalid time, incompatible
net flux and cancellation publish no partial field. These are fixed-geometry,
current-free static snapshots; histories alone do not imply induction, heat,
moving media or an infinite-domain certificate. The G1/G3/G4 regressions in
`tests/scenario_magnetic.rs` exercise IR round trips, analytical B/H and energy,
normal-sign reversal, native compatibility refusal and cancellation.

`vector_potential_boundary_at` and `solve_vector_potential3` also connect
canonical A histories to the original current-carrying Whitney edge owner.
Explicit bindings distinguish uniform patch values from authored line-average
vectors on each canonical edge. The adapter rotates polar vectors once and
forms A_average dot (x_b-x_a) in webers. It never invents a spatial profile:
uniform vectors have zero face circulation, while independent edge averages
can impose nonzero magnetic flux. Original shared-edge compatibility, current
continuity, topology, gauge and energy checks remain native.

Fixed geometry uses WORLD Cartesian coordinates. Invalid metadata, duplicate
or interior faces, nonfinite edge arithmetic and cancellation refuse. These
boundaries can also feed the original induction owner. Focused native
regressions cover rotated pure gauges, exact nonzero uniform B/H and energy
from affine line averages, shared-edge conflict and interior-face refusal.

`scenario_magnetic::induction::induction_boundary3` resolves canonical
histories from the original accepted A-phi clock to an absolute endpoint;
`step_induction3` reuses that admission for the native field step. Authored
start traces must match retained edges within the
explicit positive Wb allowance; the native candidate owns full fixed support,
charge/current compatibility and magnetic-boundary reaction work. Fixed-frame,
nonperiodic affine segments and held steps are admitted. Interior knots,
clamp transitions, endpoint jumps and dynamic frame ancestry refuse rather
than being smeared into another physical forcing law.

`step_rotor3` feeds that same admitted boundary into the original finite-inertia
rotor closure. Source/frame/trace errors and native field/mechanical errors
remain separate structured variants. The public boundary resolver also feeds
the existing joint finite-temperature and complete heat-provider methods, so
thermal drivers need not publish a standalone field or rotor step first.
The resolver itself advances no history; native receivers still admit complete
fixed support, gauge and shared-edge compatibility. A canonical rotating-field
regression compares original cochains, speed, Joule heat and boundary supply
with the direct native rotor candidate.

The adapter returns the original unpublished induction candidate and its
Joule heat, Lorentz work, magnetic storage and boundary-work account. It adds
no heat integrator or inferred thermal receiver. The original tree gauge now
preserves exactly the edges whose endpoints share one fixed trace component,
avoiding cancellation from adding and subtracting an identical gauge shift.
New regression definitions cover the analytical one-tetrahedron Faraday field
and 1/160 J Ohmic heat, exact native replay, split affine steps, stale histories,
unsupported intervals, selected frame ancestry and late atomic cancellation.
Pinned formatting/source review passed; native execution remains pending
because existing asupersync compilation exceeded the environment memory limit.

## Source-resolved speed-dependent Hertz friction and heat

`friction_heat::hertz::rate_dependent` connects the original full local-speed
curve to finite-body recoil and the original finite-cell thermal owner.
`propose_rigid_hertz_rate_heat` refreshes every typed curve parameter from the
same pinned ordered-interface claims at the trial common-face temperature.
The original Hertz pressure integral evaluates that complete `FrictionLaw` at
each local spatial/time speed. Prepared native response and the existing
Newton/finite-cell transaction are reused; only integrated physical loss enters
heat. The original exact and curved Coulomb entries retain their existing
default source and report types.

`propose_rigid_hertz_rate_heat_march` reuses the same full-versus-two-half-phase
controller and compensated whole-span kinetic/caloric/source balance. Accepted
receipts retain the complete curve parameters and source identity. Source
support, closed-normal geometry, independent native work, thermal acceptance,
work budgets and cancellation remain mandatory. A late failure leaves the
original bodies and thermal owner available for exact retry.

Four focused tests in `tests/articulated_friction_heat/hertz/rate_dependent.rs`
check analytic affine sliding/spin arrest and thermal roots, an independent
Stribeck Hertz moment from its radial beta-integral series, continuous mixed
motion and three-cell thermal feedback, and source-domain/late-failure/retry
behavior. Spatial-pressure, native-time and thermal-window errors remain
separate refinement axes. Each phase holds its solved-temperature parameters;
the march resolves their evolution. Poses, load, radius and footprint stay held.
No centre-speed fit, source reselection, drive, new friction history or thermal
store is introduced.

## Time-resolved Hertz friction and finite-cell cooling

`friction_heat::hertz::march` composes the original exact or curved native,
source and finite-cell phases over a held-geometry horizon. Full-step versus
two-half-step comparisons admit only the fine original physical states.
Individual body momentum, cell energy, net slip/spin and physical heat estimates
govern refinement. Accepted phase receipts retain original source identity and
native arrest/tail evidence. Stationary phases retain ordinary conduction
without frictional heat.

Explicit work limits include discarded thermal trials and accepted receipt
count. Any source, event, numerical, budget or cancellation refusal leaves the
complete input span untouched. The existing retained native kinetic-change
owner and compensated caloric/source sums independently gate whole-span energy;
numerical work defects are never deposited as heat.

Comparisons estimate local error rather than certified or global bounds.
Original source/event cuts, native trajectory accuracy, spatial refinement and
fixed poses, geometry, load and radius remain explicit. This composes the
existing temperature-dependent Hertz phases; it adds no alternative friction
history, heat store, drive, normal evolution or moving footprint.

Four G1/G3/G4 additions check continuous thermal fade/arrest against an
independent ODE, rest-only hot/cold transport, curved motion under continuously
changing temperature, and late whole-span refusal with exact retry. All four
and the eleven original Hertz heat tests pass in the bounded actual-source
material, native mechanics, Newton/FGMRES and finite-cell thermal build.

## Source-resolved Hertz trajectories with finite-cell heat

`friction_heat::hertz::{trajectory,transient}` compose the original exact-direction
and curved native pair trajectories with the same finite-cell/interface-temperature
Newton transaction used by the held-midpoint Hertz owner. Each window prepares
original body geometry and complete mobility once; coefficient trials restart
that captured pair and refresh the same pinned source claims at the solved face
temperature. Original physical path dissipation, including independently bounded
extinction-tail heat, enters the existing massless face once. Numerical work
residuals remain separate. Original finite cells conduct throughout stationary
dwell without additional frictional heat.

Both native bodies, thermal endpoint, source receipt and material payload leave
together after original and combined energy gates. The exact branch checks both
normal endpoint rates, which bound its linear moving phase and constant dwell.
The curved owner checks its whole-path normal-response bound. Poses, load, radius,
geometry and inertia remain held. Endpoint coefficient sampling and backward-Euler
thermal transport require window/spatial refinement; exact mechanical arrest and
numerical root convergence do not certify a continuously varying thermal law,
unique root, moving footprint or radial flash field.

Six G1/G3/G4 additions to `articulated_friction_heat/hertz` cover sliding/spin
arrest and dwell, ordinary hot/cold conduction, independent continuous thermal-fade
refinement, mixed curved recoil/heat balances, source feedback, normal
incompatibility, late failure/cancellation, work-budget refusal and exact retry.
All six plus the five original Hertz heat regressions pass in the bounded
actual-source closure; full workspace proof remains separate.

## Finite-reservoir stream exchange driving nonlinear membrane motion

`GasMembraneThermalMotion::step_discrete_gradient_with_nozzle` connects the
first moving chamber to an explicit finite, fixed-volume reservoir through the
original source-backed converging nozzle. Ideal mixtures and homogeneous
Peng--Robinson parcels retain their native throat, choking, acoustic/caloric
state, species identities, source metadata and donor/receiver inversion controls.
Different source families or unavailable cross-species interactions refuse.
Direction follows actual accepted absolute pressures on every call.

The native nozzle first transfers species and upstream total enthalpy at fixed
volumes. The resulting chamber then drives the existing discrete Gibbs
membrane trajectory, including its exact swept-volume pressure work, native
derivative, endpoint physical pressure and all support rows. An opposing sealed
gas retains its own heat, volume and work account. Prescribed heat is external
energy only; transferred enthalpy, flow work and thermalized stream energy are
not deposited again. No nominal gas constant or nonideal pressure correction is
introduced.

The reported whole-interval chamber storage includes stream and mechanical
stages. Whole-pair energy combines the original reference-safe transfer residual
with the independently gated native mechanical integration defect; species
formation datums cancel inside their source owner before rounded extensive
reports are formed. Endpoint mass is checked independently. The reservoir and
both chamber states, mechanics, time and accepted history publish only after
all source, work, energy and final cancellation checks. Private trial motion
starts without copied old history, which moves once at final publication. A
closed nozzle gives the unchanged native Gibbs state and history.

This remains a first-order transfer-then-motion split over one physical interval.
Donor-fraction and source pressure-crossing bounds refuse rather than clip mass.
Physical-time refinement is required: a small conservation residual is not a
splitting-error estimate. The nozzle opening belongs to the rigid enclosure;
resolved pipe/jet momentum, membrane impingement and reaction-force receipts
are outside this homogeneous chamber model. No open-system entropy/passivity,
spatial gas equilibration or phase-transition claim is added.

Five focused native regression definitions cover choked ideal-mixture inflow
and continued history, discharge with an independently heated sealed opposing
gas, exact PR source/nozzle/Gibbs parity, closed-nozzle identity for both source
families, and source/native-budget/whole-pair/final-poll rollback with exact retry.
They include independent finite-fixture energy checks and endpoint pressure
reactions. All five changed Rust files parse with Rust 2024 rustfmt; the two new
implementation/test modules were formatted. The complete coupled tests were
not compiled or executed because the available native dependency closure still
lacks the original material/solid/time stack.

## Finite-body thermal contact in nonlinear membrane pressure motion

`GasMembraneThermalMotion::thermal_discrete_gradient_candidate` solves one
reciprocal heat transfer between the first original gas and an explicit finite
enthalpy/phase body over the same nonlinear mechanical interval. Each trial
delegates the existing native Gibbs pressure/work trajectory and original
enthalpy entropy chart. The shared bounded heat-feedback solver closes
`Q = dt G (theta_body - theta_gas)`; no gas law, time integrator, scalar
root method or thermal storage model is replaced.

`MembraneGasSample` retains the actual ideal or homogeneous Peng--Robinson
thermodynamic mean, including its original source/transport state. The callback
also receives the actual body enthalpy midpoint, preserving latent phase
fractions. The thermal force uses the body's original `Delta h / Delta s`,
which generally differs from midpoint temperature. Work pressure remains the
native discrete Gibbs effort and is distinct from the sample's EOS pressure.
Only this new path copies the first gas's final source sample; existing
midpoint/Gibbs calls and intermediate Newton caches retain their prior behavior.

One signed Q is debited from the body and credited to the first gas. Prescribed
first-gas, body and opposing-gas heat remain separate. A sealed opposing chamber
keeps its original volume/work/heat law. Reports independently retain all
mechanical and caloric storage, pressure work, source residual, body energy,
original gas/body entropy and evaluated nonnegative Fourier-contact production.
The entropy allowance derives from native gas, body and heat-law error budgets.
Structural integration error stays in the existing mechanical-energy account;
it is never converted into heat or used to correct the gas Gibbs balance.

The joint candidate contains motion, original finite-body proposal and balances.
All root trials use one private scratch trajectory; complete accepted history
is copied only for the final candidate. Source, body, interface, work, energy,
entropy or cancellation refusal preserves both gas states, body input and
accepted mechanical history. Conductance callbacks must remain repeatable and
defer irreversible side effects until joint publication.

The body retains fixed mass/geometry and negligible-pV-work enthalpy storage.
It is an explicitly supplied thermal inventory; membrane calorics, thermal
expansion, spatial equilibration and moving thermal geometry are not inferred.
The original fixed-perimeter structural and homogeneous-gas limits apply.
Positive evaluated contact entropy does not establish temporal or constitutive
accuracy; physical-step refinement remains necessary.

Four focused native regression definitions cover hot/cold continuing motion,
independent variable-Cv/finite-body balances and insulating parity, latent
melting/freezing with a heated opposing chamber, full PR source-sample identity
and work-pressure separation, and late source/cancellation/work refusal with
exact retry. All changed Rust sources parse; the implementation and new test
modules were formatted. The full coupled tests were not compiled or executed because
the available dependency closure lacks the original material/solid/time stack.

## Discrete Gibbs pressure in native nonlinear membrane motion

`GasMembraneThermalMotion::step_discrete_gradient` connects the original ideal
and homogeneous Peng--Robinson discrete Gibbs heat/volume laws to the existing
matrix-free nonlinear membrane trajectory. The source's mechanical work pressure
multiplies the exact path-mean volume gradient; its original fixed-heat pressure
derivative and the same endpoint geometry derivative supply Newton actions.
Both first and opposing sealed chambers use their actual accepted inventories,
independent source controls and native calorics. No gas law, entropy correction,
pressure fit or replacement structural integrator is introduced.

`GasMembraneDiscreteGradientStep` retains the original
`GasDiscreteGradientStep` receipts, including independent energy, Gibbs and
entropy-transfer residuals. Work pressure is distinct from the source sampling
EOS pressure; for PR these generally differ. The reported sampling temperature
is the source's entropy-conjugate mean, while accepted source states and physical
endpoint pressure/residual remain endpoint quantities. At zero supplied heat,
each gas's entropy change is bounded by its original source allowance. Prescribed
heat contributes gas-side `Q/theta`; this is not a complete external bath/body
entropy account or a heat-transfer law.

The same native owner may continue between midpoint and Gibbs steps without
resetting mechanics, gas storage, clock or history. The original midpoint API,
report type and midpoint field names remain unchanged. Private source-typed
proposals share the original geometry, force/tangent, work gates and final
publication; no additional trajectory storage or numerical solve is introduced.
The finite-strain material force still uses midpoint coordinates. Its integration
error stays in the enclosing energy defect, separate from gas Gibbs closure and
pressure-transfer roundoff, and remains subject to time refinement.

Four focused G1/G4 regression definitions cover original mixed ideal/PR
receipts and curved force derivatives, signed heat and independent variable-Cv
entropy/storage balances in two chambers, PR work-pressure impulse versus true
endpoint readback, and source/work/energy/late-cancellation rollback with exact
retry. The changed Rust sources pass parsing/formatting; the full coupled tests
were not compiled or executed in this dependency-limited environment. These
definitions do not establish full coupled validation or continuum accuracy.

## Work-consistent gas pressure on nonlinear heated membranes

The existing `GasMembraneThermalMotion` couples each original ideal or
homogeneous Peng--Robinson parcel's midpoint heat/volume pressure to the exact
path mean of the native membrane volume gradient. The shared `fs-solid`
pressure owner integrates the changing oriented triangle areas analytically.
Its nodal force dotted with the accepted displacement transfers
`p_mid*(V_new-V_old)`; a curved volume chart no longer introduces the former
midpoint-area pressure-work defect.

Native matrix-free Newton actions differentiate this same force. Both the
endpoint-volume derivative of the source pressure and the endpoint derivative
of the mean geometry include the `q1=2*q_stage-q_old` chain rule. Original
gas calorics, pressure derivatives, heating, opposing chamber ownership and
source admission are unchanged. Physical endpoint residuals retain the true
endpoint area and pressure. The two retained geometry-force vectors replace
the former stage/endpoint gradients without adding a trajectory-sized cache.

Pressure-transfer roundoff and the full mechanical/caloric energy defect keep
their independent existing gates. No residual is subtracted, deposited as heat,
or used to rescale momentum. Material forces still use midpoint geometry, so
nonlinear structural integration error and time refinement remain necessary.
This is a spatially uniform gas with fixed perimeter and stationary supports;
it adds no fluid inertia/waves, moving enclosure, phase change or thermal law.

Existing coupled regression definitions retain original ideal/PR caloric and
tangent comparisons, verify two-gas pressure closure while mechanical error
refines, and preserve source/energy/cancellation rollback. Their full dependency
closure was not executed here. Two tests of the actual native pressure kernel
passed, including the cubic fixed-rim work identity and endpoint chain rule;
that evidence does not establish the complete coupled integration.

## Ambient cooling around modal dry contact

The finite-cell modal dry-contact transaction applies each original Robin
boundary over native exposed area minus the actual loaded footprint overlap.
Lift restores full original exposure. Migration and recontact retain the same
finite thermal histories, and source-selected friction uses the newly solved
interface temperature in the same mechanical/thermal transaction. Original
heat-transfer coefficient, bath temperature, material and caloric objects stay
authoritative; no fitted decay or extra heat store is introduced.

Shared cells require the same declared material face and disjoint loaded
footprint intervals. Overlapping active intervals and unresolved reservoir
faces refuse. Retained-film ambient coverage remains outside this route.
Geometry prepares once per physical step. The native boundary residual and
independent external-power account use identical conductances; the existing
joint rollback applies to source, coverage and cancellation failures.

The dry-contact regressions specify analytical Robin cooling with partial
coverage, lift, recontact and bow migration, plus temperature-dependent
friction feedback and a late source refusal preserving both physical owners.

## Original fluid materials in hybrid pressure and resolved interiors

`cavity_modal::interior::hybrid::HybridInterior<'a, 'c, S = GasState>`
binds the pressure region to the same original `ThermoviscousFluid` retained
by the resolved interior. Material liquids and other admitted fluid sources
carry their actual density, adiabatic sound speed, Cp/Cv and transport through
the existing geometry-verified modal interface. The resolved mixed equations
retain their original signed thermal expansion and isothermal compressibility;
equilibrium pressure is never substituted for the fluid's acoustic stiffness.
The pressure region remains the adiabatic pressure approximation with its
explicitly declared thin-wall law. Full viscous/thermal fields and bulk viscous
stress remain in the resolved region, without an inferred pressure-volume loss.

Admission checks pressure density/sound speed and optional thin-wall transport
against the same source snapshot, and resolved density/Cp/viscosity/conductivity
against its already assembled coefficients. The default gas type, original
molecular reservoirs, full pressure and resolved fields, finite solid walls and
duct convention remain on the existing owners. Preparation uses the original
resolved Cp without resampling material properties per scattering column.
Frequency, cancellation, coordinate and physical interface refusals remain
explicit. This adds no liquid relaxation, cavitation, phase change or background
flow model; modal trace convergence and the insulated, free-tangential-traction
interface assumptions remain unchanged.

Native regression definitions exercise an original liquid material card,
signed expansion, original-owner parity, fields/work, duct conjugation
and source mismatch refusals. The executor was unavailable for this extension;
native parsing, formatting, compilation and test execution remain pending.

## Resolved interior acoustic absorption on the native thermal timeline

`porous_heating::interior::InteriorAcousticHeating` binds every original
equivalent-fluid tetrahedron to one explicit existing thermal control volume.
It evaluates the prepared operator's positive local density and compressional
loss forms from the retained coherent pressure coefficients, once at preparation.
These are physical volume-integrated cycle-mean watts. The bulk-only constructor
excludes boundary absorption, outgoing power and numerical residuals; no second
volume, density, porosity or peak-to-mean factor is applied.

The complete original-order map requires unique active destinations, matching
original integrated volumes, and the exact borrowed native caloric object.
The caller owns spatial correspondence and the uniform-temperature assumption;
equal volumes do not prove overlap. Native mass, geometry, material, caloric
identity/storage basis and body/segment bindings are checked throughout an
advance. The original source pressure and equivalent-fluid owner remain borrowed.
An arbitrary supplied pressure field is not certified as an acoustic solution.

A composable native observer adds cached watts to authored sources only for
accepted thermal steps. Contacts, reservoirs, events, caloric domains, latent
storage and timestep/conservation admission remain in `fs-thermal`.
Separate mechanism and per-cell Joule receipts use the actual accepted durations
and are already included once in the native external-heat report. Candidate
thermal and observer states publish together; cancellation or refusal leaves
their inputs unchanged. Prepared storage is bounded by an explicit cell budget
and a hard 4096-cell ceiling; thermal steps perform no acoustic/material solve.

This is one-way frozen harmonic heating. Reprepare when equilibrium/material
state or coherent drive changes; the mechanism split supplies no independent
fluid/skeleton temperature model or transient acoustic startup. Native tests in
`porous_heating/interior/tests.rs` specify independent affine-field loss
integrals, authored-event accounting, observer composition, bindings and rollback.

### Original mechanical-wall absorption on the same thermal timeline

`InteriorAcousticHeating::try_new_with_mechanical_walls` extends that bulk path
with explicit receivers for every original mechanical wall patch. Each supplied
wall state retains its solved Maxwell extensions, so tiny physical relaxation
loss is not reconstructed by subtracting rounded displacement and memory.
The original consistent trace forms supply separate positive resistance and
relaxation watts; signed pressure exchange, external traction work and numerical
balance residuals never become heat sources.

`InteriorWallHeatCell` declares the exact native caloric object and complete
one-sided acoustic area represented by its thermally lumped receiver. That area
must match the original FE trace integration, including curved/high-order
geometry, and the native thermal mass must match original surface density times
that area. Volume and wall receivers are all distinct and count against the same
bounded budget. The native total surface area, which can include both wall faces
and edges, remains independent and unchanged. The caller owns actual spatial
and material correspondence and the justification for a single wall temperature;
matching area and mass alone does not establish those physical assumptions.

Preparation evaluates each original volume/wall loss form once and retains the
original pressure operator, pressure field and complete wall-state slice. All
four irreversible mechanisms and original-order volume/wall cell energies have
separate Joule accounts on the same accepted thermal substeps. Their sum is
included once in native external heat, alongside authored sources, conduction,
reservoirs and events. Binding checks and native transactional rollback cover
both kinds of receiver. Bulk-only construction remains explicit and deposits no
wall work. This remains frozen harmonic heating; no coherent field recertification,
wall-temperature feedback or new solid thermal law is inferred.

`porous_heating/interior/tests/walls.rs` specifies actual coupled pressure/wall
states against independent triangle/Maxwell heat integrals, native energy and
event accounting, contact-driven timestep admission, split observer history,
zero-drive parity, complete area/mass/caloric binding, and late refusal rollback.

### Resolved viscous acoustic heat on native thermal cells

`porous_heating::interior::thermoviscous::ThermoviscousViscousHeating` connects
the original full mixed P2/P1/P2 acoustic loss observation to the existing
`DepositionSimulation` timeline. It borrows the exact prepared frequency,
transport/caloric state and original velocity/temperature coefficient slices.
The native quadrature evaluates every retained affine or curved tetrahedron
once. Only positive shear and bulk viscous powers become caloric sources;
conductive entropy production and thermal/molecular availability destruction
are separately accumulated observations. Signed port/source work, first-order
heat forcing, boundary conduction, structural losses and numerical residuals
are not deposited by this fluid-volume adapter.

The caller supplies one distinct original native caloric receiver per element.
Both physical volume and equilibrium fluid mass must match, and the existing
binding checks retain material, caloric object/identity, storage basis, geometry
and body/segment state. A prepared element-count query admits the complete cell
budget before quadrature. Matching measures do not prove spatial correspondence
or caloric equivalence; those remain explicit caller responsibilities.

Cached watts compose with authored thermal sources, contacts, reservoirs and
events through the native observer. Accepted durations alone enter separate
mechanism and original-order per-cell accounts. Shear plus bulk heat is included
once in native external heat; availability loss never joins that total. Native
candidate thermal and observer states publish together, preserving the input
state and work account on cancellation, binding change or caloric refusal.
Continuation uses the accepted observer without another acoustic evaluation.

This is held, quasisteady cycle-mean heating, with no material refresh,
mean-flow/streaming closure, acoustic startup or moving-geometry inference.
`porous_heating/interior/tests/thermoviscous.rs` supplies independent affine
strain/divergence and thermal-polynomial oracles, native event and observer
composition, separate availability accounting, original bindings and rollback.
Formatting and bounded source review do not establish native runtime success;
native execution remains pending in the current pinned-dependency environment.

### Molecular duct viscous heat on native fixed-volume gas cells

`porous_heating::interior::duct::MolecularDuctViscousHeating` connects original
molecular cylinder and spherical-cone loss forms to `DepositionSimulation`.
Both optional forms must borrow the exact same `PreparedMolecularNetworkFields`
object; numerical equality of different prepared source banks is insufficient.
The original network, material/diffusion/wall sources, frequency, clock and
coherent mean/covariance factors remain borrowed. Each native loss bank evaluates
coherence once at construction. Independent factors add powers after coherent
interference within each factor; no default excitation replaces an absent mean
and empty factor list. No cross-frequency or finite-shutter heat is inferred.

Destinations follow cylinder cells then cone cells in original input order,
independent of CellId sorting. Every selected interval must have positive axial
length and volume, and its native receiver must have the same length, physical
volume and original gas density times volume. Zero cells, repeated destinations
and overlapping intervals on the same original edge/authored segment refuse.
Cone slice identities do not permit overlapping physical intervals to be counted
twice. A selected subset of the graph is supported; uncovered regions are not
silently credited. Absolute world placement/orientation remain caller-owned,
and the native exposed convection area is never replaced by acoustic wall area.

The same original `CellBinding` retains geometry, mass, material, body/segment,
exact caloric object, storage basis and caloric identity. Construction requires
source equilibrium temperature and an identified InternalEnergy law. The caller
must supply that gas's original Cv, datum and caloric correspondence; matching
rho/T or a supplied identity is not proof of shared material provenance. An
enthalpy/Cp law cannot stand in for fixed-volume internal-energy storage.
No caloric law or gas equation of state is fabricated by this adapter.

Only native `viscous_power_w` enters thermal sources, with no second volume,
half or absorption factor. Fourier and molecular availability destruction and
mechanical-wall work are retained as separate accepted-duration observations;
none is deposited as gas heat. The source map is prepared once. Native thermal
sources, contact conduction, reservoirs, events, material domains and timestep
admission retain their existing owner. These same viscous watts must not also
be authored as separate thermal sources.

`DuctHeatingObserver` shares the immutable source and bindings while cloning
only bounded energy accounts. Sources, events and accepted steps revalidate
original bindings and the caller's finite positive temperature-excursion limit
from each frozen gas state. Each accepted duration is charged once; there is
no final-duration reintegration. Native candidate state and observer history
publish together. Cancellation, material/binding or excursion refusal preserves
the supplied accepted state and work account. Continuation uses the accepted
observer and requires no acoustic solve, material refresh or loss quadrature.
The convenience `advance_to` starts a fresh account for that requested advance.

This is a held cycle-mean, frozen-material, fixed-volume coupling. Temperature
excursion is an admission limit, not a certificate for unchanged material/phase
or a model of pressure work, acoustic startup, streaming, thermal expansion or
feedback. Magnus-profile interiors still require their original spatial law.
The three focused native tests in `porous_heating/interior/duct/tests.rs` cover
mixed coherent cylinder/cone heat with original caloric/event accounting,
native contact/reservoir transport and split/composed observer history with
late-cancellation rollback, and source/caloric identity, geometry/mass/overlap
and strict-temperature refusal/retry. They use synthetic original-source
fixtures, not empirical validation. Original Rust formatting/parser and
whitespace checks passed. Native tests were not executed; source review and
formatting do not establish compile or runtime success.

### Accepted-temperature feedback through the full interior pressure field

`porous_heating::interior::feedback::TemperatureInteriorHeating` connects the
original native caloric temperatures to the existing `PorousTemperatureFluid`
gas/liquid/molecular sources and complete equivalent-fluid FE owner. Each original
tetrahedron retains its own fluid state and, where declared, the original fixed
`PorousFluidLaw` response. Spatially distinct temperatures are never averaged.
Initial volume-fluid sampling occurs after native FE work/storage admission;
later refreshes reevaluate only changed cells before one global coherent solve.
Original pressure, all distributed mechanical coordinates and their separate
volume/wall physical loss forms remain in that solve.

The native thermal step consumes powers from its accepted START state. Those
watts times the actual admitted duration enter the shared four-mechanism account
before endpoint material/field refresh. The returned acoustic owner is already
evaluated at the returned thermal temperatures. Unchanged temperatures, including
latent plateaus, retain the same immutable prepared factor and heat snapshot.
Cloning shares those immutable operators, while candidate work and energy
accounts remain independent. Source, field, binding, domain or cancellation
failure publishes neither a partial thermal state nor a partial acoustic cache.

This is explicitly quasisteady feedback at prescribed source pressure and
composition. The native thermal law must justify thermally lumped specimen
storage and its common fluid/skeleton temperature. Thermal mass, geometry,
caloric object and storage basis remain unchanged; mass is never replaced by
rho(T)*volume. The closure does not evolve a sealed gas inventory, background
flow, pressure or separate constituent temperatures. An ordinary free-fluid entry
uses its refreshed lossless adiabatic EOS; viscosity alone is not an invented bulk
attenuation. Original molecular sources reevaluate the complete reservoir law
and explicit composition at the new gas state, using `effective_relaxing` in
pores and the native free molecular response outside them. They never rescale
old reservoir times or add molecular attenuation to an equilibrium pore law.
The geometry-free `PorousTemperatureFluid::volume_response` hook consumes these
same cached source states. A custom source without this explicit constitutive
hook refuses, preserving dispersive behavior instead of falling back silently
to its equilibrium `ThermoviscousFluid` tuple. The existing coherent-layer hook
retains its own unchanged native owner.

The owned frequency, coherent velocity drive and Robin data remain fixed.
`try_new` also holds its explicitly authored mechanical coefficients.
`try_new_with_material_walls` instead binds each original
`CylindricalMaterialBoundary`, material card and temperature coordinate to its
complete native wall receiver. It evaluates the same source at that receiver's
accepted temperature, retaining property claims, axes, restraint, thickness and
applicability. Original Young-viscosity, source-pair Maxwell and admitted
principal-aligned laws use `LocalMaterialWall::at_temperature_pinned`; unsupported
constitutive continuations, changed density/mass, phase and expansion refuse.
The original cylindrical faceting and harmonic reduction are shared with
`MaterialMechanicalHarmonic`, rather than inferred from a replacement impedance.

Wall-only temperature changes refresh the complete field without reevaluating
unchanged fluid sources. Simultaneous fluid/wall changes share one preparation.
Unchanged wall states share their immutable full material receipts; accepted
wall temperatures, material states, mechanical coefficients, field and heat
advance together. The source constructor requires empty raw mechanical inputs
and never silently replaces an authored wall. `problem()` exposes the current
accepted mechanical coefficients, while `wall_bindings()` retains the original
sources and `wall_materials()` exposes their accepted material snapshots.

Thermal timestep refinement must resolve feedback independently of conduction
CFL; acoustic startup, parametric stored-energy work and moving geometry remain
outside this cycle-mean approximation, as does transient wall-memory transfer
between harmonic solutions. The accepted factor remains live while
a new preparation is built, adding its retained storage to the native
per-preparation memory allowance. Native substep and FE budgets remain binding.

`porous_heating/interior/tests/feedback.rs` specifies original-source field
parity at distinct temperatures, accepted-start joules with native contacts and
events, split continuation, latent/one-cell cache reuse, source-domain and
cancellation rollback, and temporal refinement against independently
orchestrated original material/FE/frozen-heat owners. Its molecular child compares
changed-temperature free-gas and porous responses with the original joint
reservoir owners and retained source composition. The wall child exercises
source-backed temperature-dependent walls, independent frozen-source heat and
field reconstruction, changed-wall work/cache behavior, and refusal rollback.

#### Coherent physical excitation during temperature feedback

`InteriorTemperatureAcoustics` retains explicit `boundary_drives` and
`traction_drives` alongside its existing linear velocity traces. An
`InteriorTemperatureBoundaryDrive` supplies the original dimensionless
six-node quadratic trace, explicit source quadrature and complex peak inward
velocity [m/s]. An `InteriorTemperatureTractionDrive` supplies the original
global hierarchical `WallTractionProfile` and complex peak outward pressure
[Pa]. All drives share the native exp(+i omega t) convention and preserve
their signed shape and relative phase. Traction drives wall momentum;
the complete pressure, displacement, velocity and Maxwell memory remain solved.

Each changed fluid/wall temperature creates one native equivalent-fluid
preparation. Its original quadratic-source and wall-traction owners prepare
against that exact factor, then `solve_with_sources` combines every drive
before one coherent full-field solve. An old prepared source is never reused
with new constitutive coefficients. The authored shapes and amplitudes stay
fixed; unchanged temperatures reuse the existing immutable field/heat snapshot
without source projection or material evaluation. There is no factorization
per drive, substitute source equation or prescribed mechanical velocity.

Existing source-count limits and finite amplitudes are admitted before material
evaluation. Native source quadrature, original wall coordinates, boundary
support, pressure-release conflicts and cancellation retain their existing
refusals. Fluid and wall material refresh, physical fields, thermal cells and
heat receipts still publish through the original common transaction.

The field retains signed boundary-source and applied-traction work separately.
The thermal timeline receives only the existing positive volume-loss and
wall-resistance/relaxation forms; source work, internal pressure exchange,
Robin loss and numerical residuals are not added as another heat input.
Accepted-start power and actual native substep duration determine deposited
joules. These remain frozen harmonic source drives within quasisteady thermal
feedback, without acoustic startup, external-actuator dynamics, transient
mechanical-memory transfer or parametric stored-energy work.

`volume_drives` extends the same transaction with original distributed
continuity injection. Each `InteriorTemperatureVolumeDrive` retains an immutable
thread-safe spatial callback receiving the original tetrahedron index and
physical position, an explicit native quadrature order, and a complex peak
volume-flow amplitude [m^3/s]. Its signed profile [m^-3] is neither normalized
nor multiplied by an inferred density, porosity or heat-release coefficient.
Clone shares the original callback; Debug never invokes it. Profile values must
remain fixed for the same arguments throughout the owner's lifetime.

Every changed-temperature factor prepares its own native volume source before
the coherent boundary/volume/traction solve. Unchanged-temperature snapshots
perform no callback or projection. Native source caps and finite amplitudes
are admitted before material evaluation; source quadrature, nonfinite profile
values, callback failures and cancellation keep their original refusals. Source
work stays separate from the existing physical heat forms. Focused regressions
specify analytic homogeneous compression/work, plateau callback/cache reuse,
mixed original-source heat at changed temperatures, admission and rollback
during profile evaluation. Native compilation/execution remains unverified
under the current missing constellation dependency.

`porous_heating/interior/tests/feedback/sources.rs` defines focused original
field/work comparisons, physical heat accounting, material refresh and source
admission/rollback cases. Native parsing, formatting, compilation and execution
were unavailable for this extension; source review alone is not runtime proof.

## Optional local material tangents in finite-solid motion

`FiniteSolidBody::with_tangent_storage` selects ordinary midpoint Newton tangent
storage without changing the accepted physical state. `Assembled` remains the
default. `LocalHessians { maximum_tangent_bytes }` retains the original 81
constitutive coefficients per tetrahedron, applies the native prepared material
operator, and obtains its signed diagonal from the same local coefficients.
The allowance covers only that Hessian payload; reusable full direction, force
and scratch vectors, a reduced diagonal, and optional external sparse
derivatives require additional memory. Insufficient allowance refuses tangent
preparation without advancing the accepted body.

Local actions avoid global material stiffness triplets and reuse their vector
storage across Krylov products at the same displacement. Prescribed-coordinate
variations remain zero even for nonzero held displacements. Follower pressure
and gravity retain their original signed, potentially nonsymmetric derivatives;
the preconditioner receives the material plus external diagonal before its
existing scalar safeguards. Local storage is an explicit choice because its
retained memory and product cost can exceed compact assembled CSR on connected
meshes. Energy-momentum integration keeps its separate assembled algorithmic
endpoint Jacobian and original endpoint/midpoint chain factor for either choice.

Finite-solid material tangent actions forward the original directional
constitutive derivative. Thermal references map increments linearly by the
original inverse reference deformation and pull Piola increments back with the
same original volume/area factors as forces; tiny variations are never formed
by subtracting two identity-bearing deformations. Original material domains,
positive Jacobians, external-load owners, candidate acceptance, and cancellation
remain authoritative.

## Original material membranes with lossy three-dimensional acoustic volumes

`heterogeneous::mechanical::MaterialMechanicalHarmonic` retains the original
cylindrical wall source and each free-fluid/porous source alongside one complete
pressure/wall solve. The shared geometry admission checks the actual faceted
cylinder, outward orientation and original membrane thickness/radius domain.
Every wall's original frequency band must admit the requested solve frequency
before any volume source is sampled. Reported wall-band intersection is separate
from volume-source applicability, which the original material owner admits.

The same `LocalMaterialWall::acoustic_cell_wall` reduction used by the transient
and duct owners supplies surface mass, equilibrium stiffness, Kelvin resistance
and every positive Maxwell arm. Real coefficients enter the native harmonic
equations directly; consistent FE trace area is applied once. The source's
negative-frequency-exponent impedance and the pressure solver's positive
exponent are not mixed. Original material snapshots, zero source arms, bands,
restraint and geometry remain available through the retained source bindings.

The generic `PorousRegionsInterior::new_with_coupled_mechanical_walls` constructor
also admits explicitly authored numerical walls. Both paths preserve complete
pressure and wall coordinates, coherent applied traction, and separate boundary,
volume and mechanical heat accounts. Original local membrane and rigid porous
frame assumptions remain binding; no shell bending, moving pore skeleton,
thermal feedback or bulk dispersive stored energy is inferred.

`tests/material_mechanical_acoustics.rs` compares the coupled porous/material-wall
field with the separately evaluated original constitutive impedance, including
all pressure coefficients, wall velocity and independent fluid/wall/total work.
It also specifies original bandwidth, geometry, source and cancellation refusals.

## Original free-fluid and porous sources in resolved acoustic volumes

`cavity_modal::interior::heterogeneous::porous::PorousRegionsInterior` binds
one original source slot to every tetrahedron of the complete pressure FE
owner. Lossless free fluids supply their coherent `FluidAcoustics` state. Rigid porous
regions use the original `PorousFluidLaw::effective_fluid`, or its joint
`effective_relaxing` molecular/pore-wall response, at the exact solve frequency.
Each used region is evaluated once after mapping and FE-budget admission;
unused slots are not queried. Sources, region maps and actual frozen responses
are retained through repeated coherent boundary drives.

The native porous and interior owners share `exp(+i omega t)`. Effective rho/K
already include porosity; no additional porosity, phase conjugation, air
attenuation or surface-impedance substitution is introduced. Pressure and
superficial normal velocity meet through the common original weak form. Local
volume heat comes from the original separate density/compression loss forms;
exterior boundary dissipation remains separate. Geometry, h/p and material
quadrature require convergence; the law remains rigid-frame, small-signal,
fully saturated and within the original source's fluid/frequency domain.
No elastic Biot waves, evolving pore state or new empirical qualification is
inferred. Reprepare for changed frequency or equilibrium.

`tests/porous_interior.rs` specifies a three-rung full-volume JCAL slab against
the original analytic layer impedance, a free-fluid/porous layer with interface
and area continuity, and exact source reuse, work closure and native refusals.

### Free molecular gas in the same resolved pressure volume

`AcousticVolumeRegion::MolecularGas` retains the original `RelaxationState` and
uses the existing `PorousFluidResponse::relaxing_fluid` constructor, shared with
the native free-gas layer owner. Every admitted molecular reservoir contributes
its original dynamic heat capacity to `K(omega)=p*Cp(omega)/Cv(omega)` while the
equilibrium density remains unchanged. That material constructor converts the
original `exp(-i omega t)` calorics exactly once into the FE owner's positive
exponent. The distinct `MolecularGas` state retains the frozen frequency and
coherent gas response, with exactly zero static resistivity and no skeleton.

This is molecular compression loss and dispersion, including uniform breathing
fields; no propagation attenuation is fitted into a modulus. Stokes–Kirchhoff
attenuation, bulk-viscosity fits, Fourier diffusion, shear and wall layers are
not added. Source evaluation/domain admission and unused-region behavior stay
with the existing owners. The same complete pressure solve, material-interface
continuity, distributed mechanical-wall states and positive element heat forms
apply. Volume heat remains separate from physical wall dissipation, and no bulk
dispersive stored energy or thermal feedback is inferred.

`tests/porous_interior/molecular.rs` specifies an independent complex reservoir
sum and constant-pressure volume heat, the no-reservoir limit, a conforming
molecular-gas/relaxing-porous mesh ladder against the native exact layer solution,
and shared molecular/mechanical-wall work with retained sources and cancellation.

## Free-free cooling in reacting magnetized gas heating

`MagnetizedIonizedGasHeating::advance_reacting_with_optically_thin_bremsstrahlung`
adds an original electron-ion continuum source to the same held Maxwell work,
atomic reaction and thermal-exchange solve used by `advance_reacting`. Both
entry points share source binding, endpoint increment limits, cancellation and
atomic refresh of original transport/background laws. The native gas result
retains the cooling source and its separate escaped energy; continuum emission
does not create bound-free photon counts or a second Maxwell heat debit.

Both reacting entry points require coherent photon energy strictly below every
supplied original atomic photoionization threshold, including absent stages.
The prepared carrier-friction dielectric law has no direct bound-free
absorption, so an ionizing coherent drive is outside this coupling closure.
Selecting the new cooling entry explicitly assumes optically thin free-free
escape. Source domains and native time/error/work controls remain authoritative;
field evolution, radiation trapping and new frequency-resolved continuum data
are not inferred from an integrated emissivity.

## Finite coherent gas films with reacting diffuse radiation

`ionized_wave_heating::layers::LayeredIonizedGasHeating` binds every original
finite Maxwell film to its exact gas state, atomic sources, complete momentum
bank and immutable bound-charge material card. All finite films must be
contiguous gas parcels with positive volume equal to area times thickness.
Lossless films require an exactly neutral gas. The coherent drive must remain
strictly below every supplied original photoionization threshold, including
initially absent stages; this dielectric law supplies carrier friction and does
not supply direct bound-free absorption of the drive.

One native multicomponent diffuse column evolves all films together. Each
film's signed electron and heavy Maxwell powers enter once, alongside original
electron/heavy energy exchange, impact ionization, capture and reciprocal
spectral recombination transport. Emission can create the first carriers in a
different film. Optional original free-free cooling requires the caller's
explicit optically thin escape assumption and remains separate from bound-free
photon counts and interfaces. The bridge independently checks actual gas storage
plus both escaped radiation energies against original Maxwell work. Signed
Coulomb pool partitions retain a separate arithmetic roundoff allowance.

The accepted endpoint must satisfy the declared heavy/electron temperature and
per-species population increment limits. An initially empty electron pool has
no physical initial electron temperature increment; newborn carriers still
require native source support and obey the absolute population bound. Every
complete transport bank and original material card refreshes successfully before
any result is returned. Cancellation, source-domain and conservation failures
leave the original column reusable without partially published endpoints.

This is a held-field interval on uniform fixed-volume parcels. Its refreshed
transport and backgrounds seed the next Maxwell preparation; the adapter does
not continuously update coherent fields or account for a change of their stored
energy between preparations. Diffuse photons retain the native quasi-static
escaping boundaries, independently of coherent dielectric reflections.
Composition-dependent applicability of original material/ion-drag cards,
spatial mixing, conduction and flow require their own supported models.

## Reciprocal heat feedback for progressive finite ductile sheets

`solid_enthalpy::finite_ductile::sheet::FiniteDuctileThermalSheet` binds each
original P1/Q1 element to one exact source-backed caloric cell while retaining
all native quadrature histories and volumes. Cell and material clocks, source
identity, original volume/mass and temperature domains are checked. Prior
plastic/damage loss cannot be adopted without its original caloric/nonthermal
history. An optional moving sheet uses the same source reference densities and
original consistent membrane mass.

The existing bounded Picard/IQN-ILS feedback driver owns the thermal timeline,
transport, events, observer replay and whole-window energy gates. Each proposal
starts from the same accepted sheet, thermal state and observer; its actual
duration drives rate flow and initiation, and its endpoint cell temperature is
applied to every original point in that element. Only the complete converged
mechanical/caloric/observer candidate escapes. Failed coupling trials do not
accumulate extra plasticity, damage, heat or observer endpoints.

Original point receipts and reference volumes integrate explicit fractions of
native plastic dissipation and damage release. Complementary fractions remain
in separate irreversible nonthermal inventories. These fractions apply to
already-accounted physical losses, not total plastic work; energetic hardening
is already retained mechanically. Signed material endpoint defects and numerical
kinetic loss never become caloric heat. Source mass is conserved through damage
and thickness changes.

Temperature affects the supplied supported initiation locus. This adapter does
not invent temperature laws for elastic moduli, yield/hardening, viscosity,
expansion or post-initiation work. Thermal cells/conductances retain reference
geometry; through-thickness thermal gradients, deformed Fourier transport,
melting, contact and fragment evolution require their own models.

Three focused source tests cover an independent biaxial heat/onset/thickness
reference on P1 and Q1 with both Picard and IQN; warped per-point quadrature,
different cell thicknesses/heat fractions, rate-dependent loading and held
relaxation with retained motion/nonthermal inventory; and source, clock,
prior-loss and failed-feedback rollback. Both new files passed pinned rustfmt
parsing/formatting through stdout without filesystem writes. Bounded API/static
review passed. Runtime and type-check execution remain unverified because the
shared workspace exhausted storage and memory; no numerical test result is
claimed.

## Complete transient Rayleigh fields in a molecularly relaxing gas

`baffled_pressure::molecular_field::PreparedMolecularBaffledField` maps each
original surface quadrature point and signed input basis to the existing
material Green-function owner with Q=2*A*sum(phi*v) and Qdot=2*A*sum(phi*a).
The rigid-baffle factor is applied once. Every original path retains its own
position, direction and calorically frozen arrival front. Original reservoir
absorption and dispersion enter pressure and the full spatial velocity gradient;
reactive velocity is never reconstructed from pressure/impedance.

One caller-owned consistent velocity/acceleration history supplies every query
on the original DD physical clock. Initially resting input coordinates and
unexcited medium are explicit; missing onset impulses or nonfinite motion
refuse. A finite acceleration jump is admitted without a smooth-onset
second-order claim. The finite-prefix record includes coherent pressure,
Cartesian velocity and signed instantaneous intensity. Ideal first-order
microphone orientation may change on that record, with the same gas's
equilibrium rho*c retained as an explicit calibration reference.

The adapter prepares one original BDF2 convolution-quadrature cache, adds DD
surface projection, and includes retained/temporary mapping arrays and
projection work in explicit budgets. This is offline prescribed-motion
observation in a stationary flat infinite rigid baffle and homogeneous relaxing
gas. Classical viscothermal loss, wall layers, finite-baffle scattering,
motion and reciprocal loading are not supplied by this path. Source, spatial,
time and contour refinement remain independent.

Five regressions cover the independent lossless retarded Rayleigh field,
multiple physical fronts, molecular pressure cancellation with surviving
velocity, an independent single-reservoir Bromwich inversion with timestep
refinement, quiet modal residuals and transactional admission. Pinned rustfmt
and bounded source/API/physics review passed. Native Rust tests remain
unexecuted because the constrained real-source build stopped in asupersync
before reaching the affected crates; numerical oracle checks are not a Rust
test pass.

## Bandlimited observation of resolved physical sources

`NsfDgPerformance` and `NsfThermoviscousDgPerformance` opt into the shared causal
fractional-rate observer with `with_bandlimited_observation(policy)`, attached
after all reservoir surfaces and material-gas receivers. The original inlet
pressure, terminal flows and every selected reservoir tile form one channel
bank before Rayleigh differentiation, propagation and superposition. Distinct
source identities, signs and physical gains remain separate. The filter's
actual band, work/storage dimensions and common numerical latency are exposed
by `bandlimited_observation()`; zero prehistory and that delay remain audible.
No future source intervals, physical tail or delay removal are introduced.

Every output still advances the original physical trajectory through exactly
`source_intervals_for` coverage, including the final partial interval. Controls,
mechanics, energy accounts and callback-independent source order are unchanged.
Legacy constructors retain their original exact or conservative-overlap path.
Cold transfer into a complete acoustic field moves the original unadvanced
filter only with identical physical dt, output window and channel basis; later
surface/gas attachments refuse. Cancellation and each owner's existing accepted
physical-prefix/poison semantics remain in force. This finite observation FIR
does not undo aliasing already present in the source or qualify source dynamics,
spatial radiation resolution or empirical instrument fidelity.

## Modal convergence with conservative and cancelling fields

The modal vibroacoustic convergence APIs now require both the original active
input-power comparison and a normalization-invariant reference energy-norm
comparison of the complex displacement/pressure fields. The latter retains
conservative resonances, source-work cancellation, omitted coordinates and
phase error. `ConvergedResponse` exposes `input_power_delta` and
`energy_norm_delta` separately; their maximum is `truncation_delta` and
controls refusal. Existing callers can therefore receive a truncation refusal
for a field that previously passed on power alone.

This is observed sensitivity to one basis reduction, not a bound on pointwise
pressure or modes absent from both bases. The fixed reference metric includes
structural kinetic/strain and fluid compression/kinetic energy; additional
storage owned by dispersive materials, boundary loads or exterior radiation is
not claimed by this metric. The comparison uses one allocation-free linear pass
over the existing solve results and adds no dense solves. Deterministic hypot
accumulation avoids squaring quiet modal coefficients before weighting.

Five new regressions cover omitted lossless structural/acoustic resonances,
a representable quiet source whose pressure squares underflow, coherent
structural/acoustic work cancellation, an independently integrated physical
field with signed cavity-basis rescaling, and exactly zero fields. The existing
power-specific regressions now inspect `input_power_delta`. These Rust tests
were added but could not be executed in this session after the execution
environment became unavailable; source review and independent numerical
fixture checks do not substitute for a native test pass.

## Synchronized transient Rayleigh field and microphone orientation

`baffled_pressure::field::BaffledField` observes pressure and all three fluid
particle-velocity components on the same accepted source trace and physical
clock. It composes the original bandlimited differentiating receiver and
first-order receiver kernels, retaining both the radiating `Q'/(c R)` and
reactive `Q/R²` velocity terms on each original surface path. Pressure
cancellation does not erase a surviving vector field. The actual medium
remains bound to each frame; a prescribed volume flow's velocity is independent
of fluid density, while its acoustic pressure scales with density.

`FieldFrame::microphone_pressure_pa` applies the original ideal `FirstOrder`
pressure/velocity projection at any new front orientation without resetting
history. This is local observation at a stationary point, not moving-receiver
propagation or sensor-body mechanics. `intensity_w_m2` returns signed
instantaneous `p*u` after coherent field superposition, including reactive
energy return; it is neither pressure-squared plane-wave intensity nor an
independent sound-power or radiation-loading estimate.

The circular-outlet front door shares the original disk quadrature, physical
volume-flow normalization and spatial-band admission. Every source substep is
observed before native antialiasing. Pressure and velocity share the explicit
derivative/fractional-delay latency, including non-dyadic rate alignment.
All four previews and SI conversions precede every history commit. Warm
stepping allocates nothing; cold preparation retains seven existing bounded
FIR banks. No new propagation or source model is introduced. The original
flat infinite baffle, prescribed small motion, homogeneous inviscid fluid,
zero prehistory and one-way loading limits remain.

G1/G3 regressions in `baffled_pressure_field_tests.rs` cover independent
retarded spherical formulas at several frequencies and rate ratios, surviving
velocity at a pressure null, signed and cycle-averaged energy flux, density
scaling, outlet parity with the original consumer, and refused-trace rollback.
These new regressions were not executed during their implementation session:
the shared execution environment became unavailable after syntax checking and
source review. No successful native build, full test pass or physical
validation is claimed by their addition.

## Distributed interior sources under the full opening radiation load

`PreparedRadiatingInterior::response_with_sources` combines original distributed
volume injection, six-node inward boundary motion and prescribed port flows
under the same complete self/mutual exterior reaction. Sources are prepared
through the retained `RadiatingInterior::interior()` and its original
`HelmholtzFe`; exact geometry/multiport identity, original quadrature and
pressure coordinates, support exclusions and source caps remain binding.
All amplitudes are peak phasors on the existing FE exp(+i omega t) clock:
volume sources use m^3/s and quadratic boundary sources use inward m/s.

The original loaded pressure factor supplies the coherent affine source waves
`b = S a + e`. They enter the existing physical port closure: driven rows
receive `sqrt(R) Q + e`, while opening rows receive `-e + Zbar e`, where
`Zbar = R^-1/2 Z R^-1/2`. The same factor then reconstructs the full final
pressure field. No per-source factorization, BEM reassembly, passive-matrix
modification or equivalent piston replaces the original source/exterior owners.
Source projection callbacks run during source preparation, not warm responses.

`source_power_w` retains its existing driven-port meaning. Native boundary and
volume source work remain separate in `interior`; `total_source_power_w()`
adds the three signed physical contributions once. Radiation work and native
boundary-layer losses retain their original forms. The complete power residual
uses total source work and is never interpreted as material heat.

The explicit-tolerance source API checks the same original opening-pressure and
driven-flow residuals. Native source admission also applies at zero amplitude.
Cancellation or refusal yields no partial response and leaves the original
source projections and loaded factors usable. Empty source lists retain the
original response path.

Three focused native regression definitions compare the mixed source field with
an independent full FE/BEM block and analytic source moments, exercise pure
source emission and coherent cancellation, and check original-owner, support,
tolerance and cancellation/retry behavior. Native Rust formatting/parsing was
available; compilation and execution were blocked by the missing required
sibling workspace manifest, so these definitions are not a runtime pass.

The existing homogeneous gas, fixed geometry, identical planar opening triangles
and P0 exterior trace limits remain. Boundary drives describe interior
actuators/feeds. Exterior panels outside the declared openings remain rigid;
direct radiation from a moving structural wall requires its own reciprocal
coupling. Prescribed injection does not infer combustion, mean-flow sources,
acoustic startup or a causal realization.

## Original quadratic feeds in source-bound causal radiation

`BoundCausalRadiation::prepare_with_boundary_sources` binds original six-node
boundary motion and volume injection to its complete pressure field and the
same admitted full self/mutual BEM realization. Sources must be prepared
through `field()`. The original exact-time owner uses their cached spatial
columns without resampling faces, discarding midsides, dividing by net flow,
refitting an exterior or replacing the native propagation algorithm.

Original prescribed linear feeds precede the supplied quadratic feeds. Held
real boundary amplitudes are inward velocity [m/s], while volume amplitudes
remain m^3/s. The native exact pressure-trajectory contraction includes their
signed source work once; there is no harmonic half factor. All pressure
coordinates, local wall loss, exterior history and same-interval reciprocal
reaction remain in the original accepted step. The fitted source's finite-band
accuracy, passivity class and held-flow exchange limits are unchanged.

These extra sources act on interior feed faces. Their support cannot overlap
the original exterior openings even for zero integral or amplitude. They
cannot silently become new BEM modes or radiating wall patches. Exact source
identity, release-face exclusions, full-state budgets and cancellation retain
their native refusals. Each preparation has independent zero exterior history;
advancing one leaves the admitted source template reusable.

The existing native BEM bake regression now also exercises a signed, purely
midside feed with zero net volume flow. An independent analytic P1 projection
checks the complete transient pressure field, boundary work and exterior
exchange with simultaneous volume injection, followed by original-owner,
overlap, fresh-history and cancellation checks. Formatting and source review
completed; the required sibling dependency manifest remains absent, so the
extended regression has not been compiled or executed.

## Distributed injection drives full structural and exterior fields

`HeterogeneousHelmholtzFe::project_source_load` supplies the original signed
volume-injection and quadratic boundary-velocity columns to complete
pressure-potential multiphysics pencils. The same native projector checks
source identity, finite amplitudes, source-count limits, complete support and
pressure-release exclusions. It retains every original pressure coordinate,
including constant compression and higher-order edge/face shapes. Source
preparation uses the original affine or curved quadrature and is not repeated
by subsequent drives.

`StructuralInterior::project_sources` and the prepared structural,
material-shell and exterior-loaded shell `solve_with_sources` methods pass
these physical loads into their existing complete factors. All coupled drives
use peak `exp(-i omega t)` amplitudes; the real spatial columns need no
conversion. Amplitudes imported from an `exp(+i omega t)` drive are conjugated
by the caller. In `M_a psi_ddot + C_wall psi_dot + K_a psi + B u_dot = b`,
`p = psi_dot` and `b` has units m³/s, with no extra density or frequency factor.
Boundary injection is prescribed relative to structural motion and passive
wall flow, not a second structural displacement constraint.

The original `average_acoustic_input_power_w` already includes the complete
signed `Re(p^H b)/2` source work once. Coherent source cancellation and negative
individual work remain valid; no source work is relabeled as heat. The same
original material strain, reciprocal interface work and full mutual exterior
reaction checks remain authoritative, with material heat, wall heat and
independently observed outgoing flux kept separate. Existing raw-load solves
are unchanged. This is fixed-geometry linear harmonic forcing of the retained
fluid, structure and exterior owners; no mass-transfer law, heat-release law,
nonlinear wall motion or new acoustic reduction is inferred.

## Simultaneous original shell, interior pressure and exterior radiation

`TensorShellInterior::prepare_with_exterior` couples the original material
shell, complete interior pressure-potential FE coordinates and closed exterior
BEM through one loaded native harmonic factor. It retains all free shell
translations and rotations and the existing opposite-transpose interior work
pair. The exterior resisting impedance enters through `[I_shell; 0_pressure]`:
`M_s u_ddot + C_s u_dot + K_s u - B^T p + Z_ext u_dot = F`.
There is no dry-shell or cavity-field inverse, fitted radiation spring or
additional acoustic mode reduction. Repeated coherent force and inward-flow
cases reuse this complete factor and the same exterior source columns.

Both sides use peak `exp(-i omega t)` phasors at the original material operator's
frequency. The exterior must retain the exact same original shell object.
Cavity-outward normals must agree with the original closed source facet
orientation; paired caller vertex permutations cannot change this physical
side check. Interior and exterior media may differ. Every closed exterior panel
remains, including supported scattering faces. An empty interior interface map
provides the explicit detached-cavity limit.

`PreparedTensorShellRadiation::prepare_acoustics` prepares the existing native
exterior impedance and fields without independently factoring shell motion.
Original material strain, complete-coordinate storage/work, physical-port PSD,
exterior topology/band and cancellation admission remain authoritative. The
additional real force map has an explicit `8*n*q` byte limit. Pressure and all
particle-velocity components are observed from the jointly loaded velocity and
the unchanged native BEM columns.

Applied structural input, acoustic source input, signed interior exchange,
material heat, wall heat and exterior work remain separate. Structural and
enclosing balances debit exterior-port work exactly once; independently
integrated outgoing flux must agree with that same source-bound field before a
response is returned. This is fixed-geometry harmonic coupling across the
original zero-offset closed midsurface. Apertures, open screens, finite-thickness
acoustic skins, causal material-memory evolution and a continuum accuracy
certificate are outside this interface's scope.

## Spatial series filaments on native material and heat networks

`incandescent::filament::PreparedSeriesFilament` assembles one original
`MaterialOhmicNode` per declared thermal cell. Each node receives its own
endpoint temperature and liquid fraction and deposits its Joule power into that
cell exactly once. Signed current control fixes a common current; signed voltage
control enforces equal segment currents and the sum of original voltage drops.
Heterogeneous material cards and caloric charts remain separate.

Preparation admits complete electrical source paths, exact original caloric
objects, fixed source conditions, segment centerline lengths and L*A volumes.
The retained graph may be reused only with the same cell identifiers, material,
mass and geometry. The caller's existing prepared heat network still owns axial
conduction, contacts, radiation, storage, nonlinear solution and physical time.
This assembly adds no second thermal solver or independent radiation debit.

Three native regression groups cover heterogeneous signed-drive operation with
axial heat exchange, differing latent fractions at constant temperature, and
source/cell/geometry/work/cancellation refusal. Formatting and source/API checks
passed; integration execution remains pending missing constellation dependencies.
Geometry is fixed, and the authored verification laws are not measured filament
data. General electrical topology remains with the existing circuit owners.

## Phase-consistent incandescent cooling and adaptive histories

`PreparedIncandescentSource` now accepts the original phase thermal environment
or an already admitted native endpoint. Electrical work, stored enthalpy and
phase-dependent radiative cooling use the same accepted state and fixed geometry.
The accepted `thermal_surface` retains the native material/emissivity sample for
the renderer's existing gray-diffuse binding; observing light adds no cooling
debit. Source coordinates and fixed conditions must match across these owners.

`incandescent::history::advance_history_adaptive` compares one original implicit
step with two actual half steps. It measures absolute local error in J/kg, so
melting at constant temperature remains observable. Accepted output retains both
native half steps and only their actual heat/work accounts; no extrapolated
state is substituted. A correction handles unequal half durations on rounded
absolute clocks. Every original drive knot, feedback limit and maximum duration
remains authoritative. All native solves, including rejected/coarse trials, count
against the explicit work budget, and failure publishes no partial history.

The error target must exceed native algebraic tolerance and state roundoff.
Step doubling estimates local temporal error; it does not certify accumulated
temperature, phase or light error. Existing fixed-step histories remain available.
Focused integration regressions cover latent electrical/radiative feedback,
accepted emissivity, full-chart Biot/source-context refusal, cooling refinement,
enthalpy-datum shifts and bounded work. They are source/formatter checked but
unexecuted here because the required constellation dependencies are absent.

## Full tensor-material shells coupled to interior acoustic fields

`cavity_modal::interior::structure::shell::TensorShellInterior` couples an
original `TensorPronyShell` to every original heterogeneous acoustic FE
pressure-potential coordinate, including constant compression. The material
owner evaluates its complete tensor-Prony/Kelvin storage and loss at the
requested frequency without first factoring the bare shell. The existing
interior owner then factors the full displacement/rotation/pressure pencil;
neither a bare structural nor acoustic inverse is required.

Explicit `ConformingShellFace` pairs bind unique original planar shell facets
to coincident fluid boundary faces. The fluid volume owns the outward normal.
Original supported P1 translation traces and their exact transpose exchange
pressure work once; all original rotations retain inertia and mechanical
coupling. No independent rotation velocity is invented at the zero-offset
surface. Full high-order pressure integration and heterogeneous density and
compressibility remain with the native FE owner. An empty interface is an
explicit detached limit. Real nonnegative wall admittance is relative to the
moving wall, as in that native owner.

`new_with_refined_fluid` additionally admits `RefinedFluidFace` bindings:
multiple distinct acoustic boundary triangles may lie inside one original
shell facet. The material shell and every original mechanical coordinate stay
unchanged while the acoustic mesh is refined independently. Original P1 shell
shapes are evaluated at the fluid vertices using scaled compensated geometry;
the existing pressure quadrature constructs the same `B` and `-B^T` used for
traction and displaced volume. Fluid and shell vertex orders are independent.
An explicit meter-valued geometry allowance bounds plane separation and
extrapolation beyond each shell edge; admitted shapes are never clipped.
The allowance does not model a physical offset skin. Each fluid face is owned
once, so duplicate bindings cannot add its work twice. Cross-facet intersections
and arbitrary nonmatching overlap require a separate intersection trace.

Both mechanical force/torque and prescribed inward flow use peak
`exp(-i omega t)` phasors and may act coherently. Repeated load cases reuse one
complete factor. A different frequency requires original material evaluation
again; this frequency-bound API supplies no frozen time-domain material or
frequency derivative. Every accepted response passes original full-cycle
layer strain admission and separate structural, acoustic and enclosing power
balances. Signed interface exchange is retained, while original material heat
and resistive wall heat remain distinct. No fitted cavity spring, modal loss,
or extra radiation debit is added.

Original material and complete-coordinate factor bounds apply alongside
explicit interface-face, trace-storage and trace-work caps. Cancellation and
invalid geometry/material/load refusals publish no partial accepted field.
The current adapter requires conforming or nested planar faces and fixed
small-signal geometry; arbitrary overlap/curved traces, offset skins,
tangential viscous fluid
traction, exterior radiation and causal dispersive fluids require their own
physical owners. Spatial FE convergence remains a caller obligation.

`tests/tensor_shell_interior.rs` contains independent P1-times-P2 traction and
full native material-response checks, analytical cavity-compression and
detached limits, coherent dual-drive/signed-work checks with independent wall
heat integration, and geometry/band/strain/budget/cancellation refusals.

## Equilibrium-relative pressure in coupled modal networks

Held-force `CoupledModalSystem` frames observe vibration relative to the complete
loaded equilibrium. The static operator is `K = diag(omega^2) + sum(k b b^T)`;
connector rest lengths supply `f_rest = sum(k r b)`. For each component's
displacement-quadrature weights `w_i = Im(H_i) omega_i`, preparation uses the
existing static Schur owner to solve `K z = w`. The runtime pressure subtracts
`z^T (F + f_rest)` from that component's original quadrature observation before
applying its pressure ceiling and the combined ceiling.

Preparation checks the original adjoint residual, admits setup work and storage,
and retains compact indexed coefficients. Runtime correction uses compensated
dot products and preallocated transfer data. Real-only observations retain their
original arithmetic, and uncoupled observations use the existing equilibrium-aware
modal owner. Unobserved stiffness islands are excluded from the static solve;
their free motion and applied forces remain part of the physical propagation.
A required elastic/free island with unresolved static support refuses observation
preparation through the existing support-rank admission.

The pressure map is installed with the original prepared operator and rebuilt
with coefficient, material and radiation changes. A refused replacement preserves
the accepted state and previous map. Mechanical coordinates, spring storage,
external work, heat and reaction-force accounting retain their existing owners.
Equal-endpoint force ramps and constant zero-frequency tones use the held-force
observation; genuinely varying ramps, tones and continuous DG boundaries retain
their distinct observation law. This remains a narrow-band modal transfer and
does not establish broadband radiation accuracy.

The `pressure_observation` module in `tests/modal_coupling.rs` contains independent
two-mode analytic amplitude, static-silence, signed load/release, parameter-update,
pressure-limit, mechanical-invariance and free-island regressions for ordinary
and prepared propagation.

## Incident sound drives original material shells in relaxing gas

`tensor_shell_radiation::open::molecular` composes the original full-coordinate
`TensorPronyShell` with native two-sided molecular screen illumination. One
boundary factor serves all physical motion columns and declared incident cases;
one loaded material factor serves repeated external force/torque vectors. The
blocked pressure jump supplies `f_blocked = -G p_blocked`, and the shell solves
with `f_external + f_blocked`. All translations, rotations, tensor-Prony branches,
Kelvin loss and mutual acoustic reaction remain in their native owners.

The total fluid reaction is independently projected from the coherent blocked
and motion jumps. Accepted responses check both `external = material heat +
total exterior work` and `external + blocked = material heat + motion exterior
work`, together with native port work and componentwise traction. Total exterior
work stays signed: incoming sound can do work on the shell. Molecular gas heat
is included in exterior work and is never mislabeled as far-field radiation.
Continuous incident pressure has no jump across the sheet and is not applied
twice. Prepared receivers reuse total pressure and Euler particle-velocity rows,
add incident sound once, and reject responses from another source preparation.

The equilibrium gas, frequency/material bands, original strain limits, actual
source/receiver clearance, complex-wavenumber boundary admission, cancellation
and explicit work/storage limits remain enforced. This is fixed-geometry,
small-signal stationary harmonic acoustics with peak `exp(-i omega t)` phasors,
the same gas on both sides, and a zero-offset acoustic sheet. It does not infer
finite skins, enclosed cavities, causal dynamics or measured material accuracy.
`tests/tensor_shell_illumination.rs` covers the lossless force-only limit,
independent original P1 traction, signed sound-driven power, coherent combined
drives and source-bound receivers. The broader shell-radiation and material
coupling validation Beads remain open.

## Source-bound incandescent heating on native enthalpy

`incandescent::PreparedIncandescentSource` holds one signed current or voltage
drive on the original `MaterialOhmicNode`. Its actual resistivity/conductivity
source segments bound positive Joule feedback once per drive. Each proposed
step passes endpoint power to `LumpedEnthalpyEndpoint::solve_with_internal_power`,
retaining the original sensible/latent calorics, convection, radiation, Biot
admission, cancellation, bisection ceiling and thermal residual budgets.

The electrical path must fill the declared thermal volume. Card identity,
complete temperature support, caloric/transport fixed source conditions and
quantity conventions must agree. Exact original electrical receipts survive
the accepted endpoint. `try_new_with_caloric_curve` also follows a declared
liquid-mass-fraction electrical source along the original caloric trajectory.
Pinned source paths retain every curve knot and surface crossing and reject
interior support holes. The cached positive dP/dh bound includes electrical
changes on isothermal latent plateaus; sourced sensible spans use their original
Cp minima rather than endpoint chords. The prepared bound belongs to the exact
caloric chart. Trial fractions come from the native numeric enthalpy evaluator,
and accepted current, resistance and power use that same final phase state.
Prescribed heat and terminal work are counted once;
rendering observes temperature without adding another cooling debit.

The exposed maximum step is a sufficient feedback bound, not a temporal-error
estimate. Mass and fixed geometry remain declared inputs. Expansion and spatial
current/temperature fields need their native coupled owners and are not
inferred. The temperature-only thermal transport owner cannot freeze a
declared evolving phase coordinate: a sourced transport table carrying that
fixed coordinate refuses. Declared or explicitly phase-independent thermal
transport remains admissible. `tests/incandescent.rs` supplies authored
analytic current/voltage cases, turn-off cooling, latent storage, source-context
refusals, phase feedback at constant temperature, chart binding and cancellation;
it makes no measured-filament calibration claim.

`incandescent::history::advance_history` consumes the existing SI `TimeSignal`
constant or previous-sample hold table for current or voltage. Every original
drive knot splits the physical step exactly: the incoming drive owns its
interval endpoint, then the outgoing drive continues from the accepted enthalpy.
Each held span prepares its chart-bound electrical feedback once; each substep
uses the native implicit thermal owner. Explicit maximum duration and work
limits bound temporal resolution, solves and retained output. Additional
prescribed heat is partitioned across represented intervals as a constant power.
Cancellation, source refusal or exhausted work publishes no partial history.
Energy closure does not claim temporal convergence, automatic waveform
averaging or a reconstructed signal between held samples.

## Coherent crystal absorption on native thermal cells

`layered_wave_heating::substrate::crystal::PreparedCrystalWaveHeating` binds
the original prepared uniaxial electromagnetic field to a finite planar thermal
control column inside that same homogeneous crystal half-space. Cell area,
depths, material-card identities and fixed query context come from the original
thermal body. Each cell receives the analytic integral of the coherent ordinary
and extraordinary fields in the original dielectric loss tensor. Cross terms
survive; cell-centre sampling, independent mode-power addition and subtraction
of nearly equal boundary powers do not determine heat.

Preparation is linear in cell count with bounded work/storage and cancellation.
Accepted intervals reuse the cached watts and the existing native thermal
solver. Authored heat, capacity, conduction and boundaries retain their owners;
the additional optical energy enters exactly once. Entering, leaving and
absorbed optical energy are reported separately, and the signed numerical
residual is never redistributed as heat. Failure leaves the supplied simulation
unchanged.

This is held monochromatic, cycle-mean heating with fixed optical coefficients
and optic axis. The finite thermal column does not introduce a second optical
interface, a finite crystal-film interference model, lateral temperature
variation or temperature-dependent optical feedback. Focused native regressions
compare coherent cell heat with independent depth quadrature, verify ordinary/
extraordinary interference and subdivision invariance, advance actual cell
temperatures against an independent Fourier/caloric reference, and cover
lossless illumination, source identity, budgets and cancellation.

`PreparedCrystalIlluminationHeating` extends the same native column to an
explicitly uncorrelated bank of prepared crystal fields. Each factor retains
its own coherent ordinary/extraordinary interference; expected heat adds across
factors, including distinct frequencies and incidence directions. Physical
peak V/m amplitudes contain any square-root statistical weights. No spectral
bin width, hidden normalization or inter-factor coherence is invented.

The original electromagnetic illumination and original thermal body are
borrowed. Preparation costs O(cells × factors), with four modal-pair work terms
per factor for every cell and boundary query; accepted thermal steps reuse
one cached power per cell. An empty bank supplies darkness and preserves
authored thermal heat. Typed cancellation propagates from inner optical queries.
This is a held cycle-mean source: nonlinear temperatures are driven by expected
heat, rather than asserted to be the ensemble mean temperature. Focused native
regressions compare mixed-frequency cell heat with independent depth quadrature
and a Fourier/caloric temperature reference, rotate an unpolarized factor basis,
and exercise darkness, work admission and interruption inside factor queries.

`PreparedCoherentCrystalHeating` binds an original same-frequency coherent
illumination to a rectangular footprint and the original planar thermal cells.
Each native depth cell receives its integrated physical box absorption,
including interference between incidence directions and both crystal modes.
The footprint must have the original cell area; material-card identities,
fixed optical query context, cell IDs, depth widths and accepted thermal
clock retain their native owners. The optical volume is a control volume
inside the original homogeneous half-space, not a new interface or film.

Preparation costs O(cells × modes²), with explicit mode-pair and peak-workspace
admission and cancellation during optical work. It retains one absorbed power
per cell; accepted steps reuse those watts through the common native thermal
solve. Authored heat and the coherent optical joules enter once. No cell-centre
sampling, incoherent addition across beams, or boundary-power subtraction
determines the deposited heat.

Signed outward powers and accepted-interval energies use face order
[xmin, xmax, ymin, ymax, zmin, zmax]. Lateral faces sum across depth cells;
only the first entrance and final exit faces belong to the exterior account.
The independently evaluated residual is minus the six-face sum minus absorbed
heat, and is reported without redistribution. Crossed beams can transport
energy laterally; a two-face account is insufficient for this footprint.

This is held monochromatic cycle-mean optical heating with one temperature
per transverse thermal section. Lateral temperature variation, lateral thermal
conduction and temperature-dependent optical feedback remain outside this
owner. Regressions cover finite-footprint bright/dark fringes, an independent
local-field volume/six-face integral, original depth refinement, authored-only
darkness, native temperature evolution, work/storage limits and cancellation.
Native Rust execution is required before claiming those regressions passed.

`PreparedPartiallyCoherentCrystalHeating` extends that finite-footprint adapter
to an original `AbsorbingCrystalPartiallyCoherentIllumination`. Each factor is
a coherent angular field with all beam/eigenmode cross terms; independent
factors add expected powers and may use different physical frequencies.
Their peak field amplitudes already contain square-root statistical weights.
No incoherence is imposed inside a factor, and no hidden spectral bandwidth,
second weight or deterministic beat term is invented between factors.

Both adapters share the same private native-cell geometry, material admission,
six-face reduction and single heat-debit implementation. Work admission counts
cells times the sum of each factor's squared transmitted-mode count; the EM
owner separately bounds factor traversal, including empty factors. Storage
includes the public owner and cell watts plus peak sequential EM-query or
thermal-adapter scratch. Repeated accepted intervals reuse cached mean watts.

The source, exact interface, thermal material cards, fixed query axes and cell
IDs retain their original owners. Expected dielectric heat drives the native
deterministic thermal candidate; for nonlinear laws this does not assert an
ensemble-mean temperature. Face flux and numerical residual stay outside the
heat additions. The held optical state, planar depth thermal resolution and
atomic cancellation limitations of the coherent adapter remain unchanged.
Regressions vary coherence from anti-phase to incoherent to in-phase at fixed
beam power, compare analytic footprint visibility and native temperatures,
combine distinct optical frequencies, preserve all six face energies, and
exercise exact budgets, material refusal, darkness and cancellation.
Native Rust execution is required before claiming those regressions passed.

### Finite-film angular interference on native thermal cells

`PreparedCoherentCrystalSlabHeating` and
`PreparedPartiallyCoherentCrystalSlabHeating` bind original finite-film
illumination to the shared physical-box thermal-cell cache. All four internal
waves of every beam retain the original entrance/exit faces, material tensor,
optic axis and source phase. Every cell box must lie within the original film,
including empty illumination. Footprint area, depth subdivisions, material card
and fixed query axes must match the native thermal body.

Work is cells times `(4 * components)^2`, summed across independent groups.
Positive volume-integrated dielectric heat is cached in original cell order;
all six signed face powers remain a separate optical account. Accepted native
intervals reuse these watts and add them once alongside authored heat. Storage
includes retained cells and peak optical-query or interval scratch.
Refusal/cancellation preserves the initial thermal state.

Independent groups add expected heat without invented coherence or a second
statistical weight. For nonlinear thermal bodies this is not an ensemble-mean
temperature prediction. Optics/statistics stay held; lateral thermal resolution,
carrier transients and thermo-optic feedback are not inferred. Focused
regressions combine actual film standing-wave quadrature with analytic angular
fringes/covariance visibility and native temperature/energy checks, plus exact
exit-face cells, dark out-of-film refusal, material, budget and atomic
cancellation checks. Native Rust execution remains unverified for these controls.

## Full-coordinate material shells with two-sided open radiation

`tensor_shell_radiation::open::OpenTensorShellRadiation` binds the original
`TensorPronyShell` to the existing open-screen pressure-jump BEM. Every free
translation and rotation, consistent mass, tensor relaxation branch and Kelvin
loss remains in the native mechanical factor. Actual facet-local P1 normal
velocity and its exact triangle-mass transpose exchange work; the pressure jump
already represents both faces. Supported facets remain in the acoustic mesh.
Conforming planar quartering preserves source topology and the original P1
trace, without welding, capping or replacing the shell by modal radiators.

One boundary factor prepares all physical velocity columns at a frequency;
one loaded material factor then serves repeated force/torque vectors. The same
accepted velocity drives pressure on either side, the pressure difference and
the outgoing angular flux. Independent native port, pressure-jump and far-field
work checks enforce input = intrinsic material loss + radiation. No second
radiation debit, fitted damping, pressure clipping or mode removal is applied.

The scope is forced harmonic small-signal motion with peak exp(-i omega t)
phasors, zero-offset midsurface acoustics and the same homogeneous fluid on
both sides. Structural thickness still enters material stiffness and inertia.
Finite-thickness skins, cavity-backed faces, nonlinear motion, causal time
realizations, loaded-pole benchmarks and continuum/experimental accuracy are
not established by this adapter. Original source bands, strain limits, screen
topology, wavelength, singular quadrature, angular and explicit work/memory
admission remain active. The focused G0/G1/G3/G4 integration regressions are
`tests/tensor_shell_radiation_open.rs`; this slice does not close the broader
shell-radiation validation bead `frankensim-irp-v5-instrument-realism-pexfp.106`.

## Canonical electrical histories on native material fields

`scenario_electric` binds the existing typed `Electrics/ElectricPotential`
and `Electrics/NormalCurrentDensity` rows to original exterior tetrahedral
faces. All rows use one explicit physical sample time and the original scalar
payload interpolation, outside-domain, frame and reference policies. Current
direction is explicit; A/m² is never reinterpreted as total amperes.

The lowered `SpatialConductionBoundary` feeds the original source-bound DC
solver or existing `SpatialConductionThermalNetwork`. Conductivity, exact P1
traces, component gauges, local Joule heat and contact heat remain with their
native owners. Overlapping/empty/non-exterior patches, incompatible metadata
and cancellation publish no boundary or field. These are fixed-geometry DC
snapshots or explicitly held electrothermal drives, not an implicit continuous
history integrator, motional field or dielectric/inductive transient claim.
The `spatial_conduction_enthalpy::scenario` regressions exercise canonical IR,
heterogeneous conduction, signed current and the original finite heat feedback.

## Source-bound enthalpy endpoints with prescribed support motion

`SolidEnthalpyState::advance_transient_with_supports_candidate` and
`advance_transient_with_support_accelerations_candidate` pass source-resolved
thermal fields into the native thermal/support Newmark candidates. The existing
specimen, original property-source, reference-inventory, mesh-frame and exact
shared-clock checks apply identically to all three transient endpoint routes.

The returned state retains support displacement, velocity, acceleration,
refreshed thermal stress and reactions together. Complete stiffness stays fixed
on the support path; reference inertia stays frozen. Signed thermal work and
support work remain distinct, with no inferred caloric feedback. Endpoints
still require explicit admission of intervening thermal events and timestep
refinement. The focused `solid_enthalpy` regression advances cold/warm/hot
source states through both endpoint conventions and checks clock/source refusal
without altering accepted thermal or mechanical state.


## Source-bound stator and cage heating on native induction machines

`circuit_enthalpy::motor::induction::InductionThermalNetwork` binds the original
seven-state squirrel-cage machine to distinct finite stator and rotor material
cells. Each proposed caloric endpoint resolves both pinned conductor sources
and returns their actual temperature/phase-dependent resistances to the native
flux equations. The same accepted four fluxes, relative angle and two shaft
momenta restart every trial. One native workspace is reused across all thermal
probes; only its checked winding-resistance source rows change.

An explicit positive winding conversion maps the whole-cell conductor
resistance into the power-invariant stator or referred rotor basis. It does not
alter receiving volume or multiply the already converted copper heat. Native
stator copper, rotor copper and bearing losses enter their declared finite
receivers exactly once. Signed electrical/shaft work and original magnetic and
kinetic storage close the combined balance; torque and heat are never fitted
from that residual. Existing material Fourier links and thermal boundaries
participate in the shared enthalpy endpoint solve.

Source identity, path geometry, winding conversions, receiver correspondence
and clocks remain pinned. Source/domain, event, budget, energy and cancellation
refusals publish neither machine nor caloric history. Magnetic storage, inertias,
pole pairs and bearing drag are fixed, so resistance feedback introduces no
magnetic parameter work. Endpoint thermal coefficients make this generally
first order. `bind_with_magnetizing_curve` retains the original native bounded,
odd-symmetric magnetizing current/linkage source through every thermal probe
and accepted interval. Reversible saturation changes current and stored energy
together; only native copper and bearing losses are deposited as heat. The
curve is in the power-invariant referred winding basis, not an unconverted B-H
dataset or a terminal inductance. Temperature-dependent magnetic curves,
hysteresis, skin effects, resolved bar networks, finite power converters and
experimental calibration are outside this reduced adapter.

`tests/induction_motor_enthalpy.rs` compares both winding temperatures and native
currents with an independent coupled two-coil/caloric solution, including
conduction and cooling; it exercises generating electrical return, finite shaft
work, bearing heat and complete rollback/retry. These authored regressions do
not establish the full induction-motor flagship or measured-machine accuracy.
The saturation regression also checks the original physical leakage/core
energy integral, the linkage law, Faraday-derived mean copper currents, both
temperature-dependent resistances, source-support refusal and unchanged inputs.
## Complete heat assemblies on current and voltage motors

Current- and finite-inductance voltage-driven material motors accept
`MotorHeatNetwork<A: HeatAssembly>`, deriving their original material cells
from the exact prepared heat owner. Complete thermal-contact, integrated-layer,
closed-radiation and composed reservoir/finite-link source evaluations
participate in the same native motor/enthalpy endpoint residual. Winding and
bridge heat enter their authored cells exactly once, including an explicitly
homogenized common receiver.

Full associated source samples and typed refusals survive. Graph electrical
ports and prepared contact/emissivity laws with unowned field dependence
refuse, and the existing corrected combined-energy gate admits the closed-radiation numerical net.
The original cached rotor/current owners, material winding claims, source
limits and joint publication retain their existing responsibilities. Original
geometry, fixed mechanical/L/Kt coefficients, source domains, midpoint-current
loss quadrature and first-order thermal limitations remain.

Focused regression definitions cover independent contact/current/voltage
balances, original-entry recovery, and typed late failure/cancellation rollback.
A cross-crate regression uses the actual finite-coolant convection assembly:
changing liquid viscosity reaches winding resistance, source-limited current,
rotor motion and both finite caloric bodies through the same residual.
Native execution remains unverified because the execution workspace became
unavailable before the targeted tests could run.

## Translating supports with omitted-mode contact compliance

Fixed-map normal support translation composes with the original omitted-mode
compliance in `MultiFrictionModalSystem`. The accepted contact closure is
`N^T q - G R - s`; each candidate uses its proposed support offset and the
same original contact potential, secant reaction, retained modal response and
joint friction solve. Translation and compliance may bind in either order.
Changing attachment maps retain their existing refusals. Bound nonlinear
storage composes with the compliance on fixed supports only: beside a
prescribed translation it refuses, in either binding order.

The native account retains prescribed support work `-R delta_s`, omitted-mode
storage `R^T G R / 2`, and the declared numerical loss
`delta_R^T G delta_R / 2` as separate terms. Numerical loss is not material
heat. Held offsets, original material-source admission, released preload and
joint publication keep their existing owners; a refused sample publishes no
motion, contact reaction, offset or clock change.

Focused G0/G1/G4/G5 regression definitions cover independent force/work/storage
relations, both binding orders, opposed supports, held offsets, zero-motion
parity and refusal/retry. Native execution remains unverified because the
execution workspace became unavailable before the targeted tests could run.
The fixed-map quasistatic omitted-mode approximation and time-resolution
limits remain unchanged.

## Complete force histories at moving acoustic receivers

`pcm_wav::observation::moving::force::MovingForceRenderer` now accepts the
original direct and reflected moving source/receiver owners through
`with_moving_receiver` and `with_moving_receiver_reflections`. The existing
accepted C1 volume/force histories, causal source clocks, receiver identity,
whole-window bandwidth admission and `PressureRenderer` output remain the
owners. Callers supply accepted histories and explicit zero force prehistory,
rather than implementing the past-force integral themselves.

For each arrival, every source/image evaluates its complete force past at that
image's CURRENT virtual microphone position. Only the immutable original
force/geometry knot partition is shared between arrivals; a spatial integral
from another microphone event is never reused. Fixed fluid-frame microphones
retain their cached-prefix fast path. Moving queries allocate no buffers and
charge every integrand evaluation, including failed attempts, to the original
finite lifetime work budget. Complete blocks and their receiver clock publish
atomically.

The original continuous source-speed bound resolves narrow history peaks and
checks separation from earlier source positions at the current microphone,
which is stronger than clearance of the current causal pair alone. Unresolved
past clearance, missing geometry, quadrature depth/work exhaustion and
cancellation refuse the complete block. Adaptive Simpson discrepancies remain
error estimates, not interval certificates or evidence of producer interpolation
accuracy. Both force and force derivative must vanish at the first strength
knot; unknown initial acoustic fields are not inferred.

History integrals, pressure and polar-vector velocity retain DD accumulation
through all sources and image transforms before total instantaneous `p*u`
and ideal first-order pickup. Receiver speed changes the sampled source clock,
not pressure gain or fluid perturbation velocity. The original fluid frame,
including admitted tangential uniform flow, is retained; pressure-work flux
excludes mean-flow advection. Focused regressions compare current-position
histories against analytic reflected arrays and independent dense quadrature,
retain weak force fields after volume-source cancellation, and check complete
past clearance, exact block replay and transactional work refusal.

## Tensile admission from original ribbon material histories

Fiber tensile admission composes sourced Kelvin–Voigt stress with each original
Maxwell arm, using actual endpoint velocities, nonlinear axial relaxed-strain
memory and exact linear bending memory. Signed per-layer coefficients combine
before conservative span envelopes. Original source arms, bending projections,
stiffnesses and relaxation times are checked when the native bank is bound.

Cold binding, restore and converged native candidates share this guard. Refusal
publishes no mechanical, history, energy, heat or clock change; prepared buffers
avoid new per-step allocations. Preload remains nonrelaxing. No rupture, stress
clipping or additional material evolution is inferred. Focused regressions
cover axial and bending history overload, retained-memory restore, signed
elastic/viscous/Maxwell cancellation and width partition invariance.

## Exact held-flow window-radiation damping

`WindowRadiationLoad::new_with_sampling` explicitly selects
`WindowRadiationSampling::ExactHeldFlow` for the existing cylinder-window
parallel RLC realization. The original constructor and omitted file selection
retain `BilinearMidpoint`. File-driven voices accept
`mouth_reaction window-cylinder exact-held-flow` or
`mouth_reaction window-cylinder-piston exact-held-flow`.

The continuous fit, static-mass subtraction and radiation geometry are unchanged.
Exact cell histories follow their continuous dissipative poles under the solved
held mouth flow. The implicit mouth equation uses continuous interval-mean
pressure, and the radiation ledger uses independently integrated resistor heat.
The reported sampled resistance and discrete band error use that same mean
pressure transfer; they do not reuse the bilinear frequency warp. Preparation
is shared immutably by candidate histories, with no per-sample allocation or
transition rebuilding. This forcing declaration resolves no variations inside
the held interval and adds no claim about unmodeled exterior geometry.

## Nonlinear electrostatic volume forces on native solids

`electrostatic::spatial_nonlinear_coulomb_volume_load` transfers the accepted
native tetrahedral nonlinear electric field into the existing solid P1 load
owner. It shares the original linear adapter's exact mesh/frame/connectivity
admission and integrates only the selected original free-volume-charge term
rho E. Each nodal contribution is Volume*rho*E/4; original cell order, force
resultant, first moment, affine-displacement work and fixed-DOF reactions are
retained without remeshing or a separate electric solve.

The source field retains its full material response and all sources of E.
Only free volume charge assigned mechanically to the selected material enters
this adapter. Electrode/interface/polarization forces and electrostriction are
separate contributions; an overlapping enclosing stress load must not count
the same charge again. This is an instantaneous load snapshot, with no implicit
charge transport or field/mechanical time integration. Cancellation and mesh
refusal leave the field reusable. Synthetic nonlinear charged-tetrahedron
regressions specify force, torque, affine work and transactional refusals.

## Applied magnetic boundary histories with finite induction heat

`SpatialMagneticThermalNetwork::step_with_boundary_candidate` and
`spatial_induction::heat_assembly_with_boundary_candidate` drive the complete
original tangential A boundary through the same retained three-dimensional
A–phi and enthalpy histories. Each temperature trial refreshes the original
conductivity/permeability sources and includes the prescribed boundary -dA/dt
in charge continuity and conductor-frame Joule heat. Signed magnetic reaction
work remains external electromagnetic supply; only physical Ohmic loss and
the existing reversible magnetic entropy exchange reach caloric cells.
Both native and physical total-energy gates include that supply once.
Accepted candidates retain the new trace for subsequent ordinary hold steps.
Geometry, boundary support and gauge stay fixed, and refused trials publish
neither field nor thermal state. The bulk and complete heat-provider paths
share the same residual and source-temperature feedback.

## Finite-temperature feedback on the three-dimensional induction rotor

`SpatialMagneticThermalNetwork::step_rotor_candidate` and its
`step_rotor_with_boundary_candidate` variant advance the original spatial
field, native finite-inertia rotor and finite enthalpy together. The
linear material/entropy residual is shared with standalone induction and
physical electrode circuits. Every temperature proposal resolves original
conductivity and magnetic sources, then restarts the same accepted field and
mechanical histories for the bounded speed closure. Source-weighted operators
are prepared once per temperature proposal and shared across its field trials.

The joint driver includes magnetic free-energy change minus original parameter
work, native kinetic-energy change, impressed supply, magnetic-boundary supply
and external shaft work.
The original magnetic-caloric owner deposits reversible field exchange and
recovers magnetic internal energy. Actual conductor-frame Joule heat enters
the original finite body once; Lorentz work cancels internally against rotor
motion. Final admission independently closes magnetic internal energy, kinetic
energy and finite caloric storage against external electrical/magnetic/shaft work and
external heat. Bearing reactions retain the native full impulse account.

All inner speed/field trials consume the original thermal work budget. Invalid
caller controls are admitted before remaining-budget restriction. Cancellation,
source refusal and exhausted solves publish no field, angle, speed or thermal
state. Complete heat-assembly entry points retain contact/layer/radiation and
finite-coolant source payloads and typed errors; closed-radiation roundoff is
not external supply, and unowned electrical contact/resistor work refuses.
The boundary-driven bulk and complete-provider variants retain the same
initial trace for every temperature/speed trial and prescribe the complete
original trace at the endpoint. Subsequent ordinary steps hold that accepted
trace. Boundary reaction work enters supply once and is never relabeled as
Joule heat or added to the internal Lorentz transfer.

The rotor remains a homogeneous isotropic cylinder/tube with uniform scalar
permeability throughout the field domain. A whole-body finite caloric cell
provides its uniform temperature; a nonuniform rotating material source refuses.
Endpoint thermal coefficients make the combined method generally first order.
Moving interfaces, transported magnetic memory, phase transitions and thermal
expansion remain outside this model. Focused source-bound regression
definitions cover heating-dependent braking and applied-field forcing through
two retained intervals, independent current/field/kinetic/enthalpy balance,
agreement with the complete heat assembly and trial-budget rollback.
Native execution remains unverified while constellation dependencies are absent.

## Native three-dimensional induction with finite rotor inertia

`spatial_induction::SpatialInductionRotor` closes the original Whitney A–phi
field against one native centered, axisymmetric revolute-z body. A bounded
implicit iteration uses the same accepted field and momentum on every trial;
its converged speed supplies `v = omega cross (x - axis)` and its original
current supplies both conductor-frame Joule heat and the full Lorentz wrench.
Three-dimensional charge continuity and end effects remain in the native field.

`step_with_boundary_candidate` also carries an applied affine tangential A
history through the same speed closure. Every trial starts at the accepted
field and prescribed trace; native Faraday charge balance and boundary
reaction work are retained. `step_with_materials_and_boundary_candidate`
shares the original endpoint source preparation across those speed trials.
Complete original boundary support, relative gauge and source identities are
unchanged. An applied rotating field can therefore induce both finite rotor
motion and conductor-frame Joule heat without an invented torque law.

The `induction_rotor` stage is shared with the original planar rotor. It calls
the original articulated forward kinematics and constrained inverse-inertia
map, retaining shaft input, full body force/torque and ideal-bearing reaction.
Actual load work, kinetic work, linear/angular impulse and combined field/rotor
energy must pass before the joint candidate publishes. Coupling exhaustion or
cancellation leaves the accepted cochains, source, angle, speed and time intact.
Geometry bindings and the mechanical model are shared across accepted steps;
a single velocity buffer serves every speed trial of an interval.

Admission requires homogeneous scalar conductivity, uniform scalar induced-
linear permeability throughout the field domain, and no impressed current in
the rotor. The original fixed tetrahedra approximate a complete circular
cylinder or coaxial tube whose analytic occupancy is rotation-invariant.
Constructor checks include connected oriented cells, closed manifold boundary,
full facet/cell radial bounds, end coverage and bounded chord error. They do
not certify arbitrary intersections or field convergence. Noncircular moving
interfaces, anisotropic rotating material, hysteretic-memory transport and
magnetization-interface loads require additional owners.

Focused G0/G3/G4 regression definitions cover transverse braking in both
directions, independently integrated P1 force/torque/Joule and magnetic energy,
shaft drive followed by coasting, applied-field startup, direct native replay,
retained boundary holds, bearing reactions and transactional refusal.
The fixture uses real material-resolved spatial fields and native inertia;
synthetic coefficients do not establish experimental machine validation.
Pinned formatting/syntax and whitespace checks pass. Native compilation and
regression execution remain unverified while the seven constellation sibling
repositories are absent.

## Reciprocal heat and finite porous compaction

`solid_enthalpy::finite_plastic::porous::FiniteThermalPorousSolid` couples retained non-isochoric plastic history to the original `DepositionSimulation` caloric and transport owner. A native finite porous tetrahedral equilibrium is solved at each proposed endpoint temperature; original yield heat plus signed reversible caloric exchange enters the thermal timeline exactly once. Every mechanical and thermal proposal starts from the same accepted endpoint. The complete thermal, structural and observer candidate is accepted together, after the existing temperature, per-cell heat, thermal-energy and combined-energy gates. The original bounded Picard/IQN-ILS driver is shared with dense finite plasticity and sheet forming.

`FiniteThermalPorousMaterialBinding` resolves the complete source-bound thermal porous family and preserves the exact original caloric chart, typed context, pinned Cp/density receipts and same material-card identity. Common caloric admission is shared with the original J2 owner: an explicitly selected fully solid, stress-free zero-pressure reference chart, its complete temperature interval, and a positive conservative lower bound on caloric plus mechanical heat capacity. The bound is a constitutive admission check, not proof that the nonlinear coupled return has a unique solution.

The original reference tetrahedral volume and source density fix thermal inventory and optional consistent inertial mass. Actual density uses `rho0/det(F)`, with retained plastic volume available separately. Caloric density metadata at hot states does not replace original mass. Shared mechanical reports retain physical yield heat, signed reversible exchange, internal/free-energy changes, coefficient storage and signed endpoint-work defects. Kinetic change and backward-Euler loss are accounted independently; neither numerical defect nor kinetic loss is deposited as heat.

Thermal cells and conductances remain a supplied reference network, retaining the original observer sources, transport, events and histories. No moving-area Fourier operator, pore-fluid pressure, phase flash, damage, rate dependence or calibration is inferred. The native dry associated law, affine thermal envelope, P1 locking and physical interval/refinement limits remain explicit.

The four `finite_porous_enthalpy` regressions specify a loaded adiabatic crush against an independent hydrostatic/Voce/caloric oracle, two-cell physical heat redistribution, retained hot inertial history through late refusal/retry, and original source/chart/capacity/inventory admission. Native execution remains a separate verification requirement.

## Native electrode circuits and coupled finite-temperature feedback

`spatial_induction::terminals::circuit::SpatialTerminalCircuitThermalNetwork`
advances native circuits driving nonlinear or hysteretic solid electrodes with
the original finite-enthalpy state. The nonlinear and hysteretic child modules
reuse their existing field thermal residuals; no additional thermal integrator
or approximate conductor law is introduced. The joint candidate retains the
native circuit, complete field/magnetic memory and full thermal source receipt.

Each endpoint trial resolves bulk scalar/tensor conductivity from its original
material. Optional complete graph-resistor bindings resolve every authored R(T)
from a separate whole-cell finite body. Accepted source identity, electrical
basis, geometry, caloric temperature and original resistance must agree, and
trial refresh cannot replace them. Numeric authored resistors retain explicit
heat destinations. Graph and field histories restart together on every trial.

Actual local bulk Joule and irreversible pinning heat enter the original field
receivers; authored graph-resistor heat enters its explicit receivers once.
Canceled numerical reference heat and numerical magnetic dissipation never
enter caloric storage. Stable graph/magnetic storage, finite thermal storage
and numerical loss close against original external/impressed supply and external
heat. The nonlinear path preserves its optional source-separable magnetic
entropy and parameter-work account. Hysteretic fits remain constant throughout
the admitted caloric interval, while their vector memories remain retained.

Complete heat-assembly adapters preserve original contact, layer, radiation and
finite-coolant payloads, typed failures and closed-radiation accounting.
Unowned electrical contact work refuses. Nested circuit/field trials consume
the shared thermal budget, and cancellation or any failed source, solve or
balance leaves all accepted histories unchanged. These adapters describe
stationary, fixed-geometry, fully solid bodies within the existing MQS scope.

## Actual-temperature feedback into finite-strain membrane relaxation

`membrane_thermal::visco::temperature::advance_candidate` feeds each original
caloric receiver's actual temperature into native Maxwell relaxation times.
`ViscoMembraneKinetics` retains the original unshifted spectra and one explicit
WLF/Arrhenius calibration per triangle. Every preparation uses the shared
`TemperatureShift::shift_hyperelastic` owner; repeated intervals never rescale
an already shifted spectrum. The native owned-motion kinetic gate pins all
energy parameters, natural references, branch ordering and solver controls.

The original thermal timeline chooses the interval. Its accepted starting
temperatures set the kinetics; native physical branch loss returns to the same
enthalpy cells. Accepted steps and thermal events synchronize endpoint times
without changing history, physical storage, momentum or interval heat. The
frozen Kelvin--Voigt, frozen Maxwell and kinetic-feedback routes share receiver
admission, source/inventory checks, physical/numerical work accounting and atomic
candidate publication. Numerical kinetic damping never heats a cell.

Both kinetic calibration support and the explicitly admitted fully-solid
receiver interval are checked, including after final heating. Cancellation,
foreign energy/reference laws, out-of-domain endpoints and late observer
refusals preserve all accepted owners. Temperature-independent free energies
make this kinetic feedback compatible with the existing caloric account; no
reversible thermal work is inferred. This is first-order lagged-temperature
coupling, with explicit timestep refinement still required. Thermal expansion,
temperature-dependent elastic moduli and reciprocal entropy exchange remain
separate owners.

`tests/visco_membrane_thermal/temperature.rs` covers native shifted/reduced-time
agreement, self-heating at three timestep rungs, heterogeneous WLF/Arrhenius
regions, event/continuation equality and late whole-call refusals.

## Retained kinetic storage through full and partial rigid liquid release

Rigid liquid release compares original native momentum/inertia storage with
the independently integrated P1 field at a fixed COM velocity reference,
retaining mass and momentum defects. Its final gate adds the native transfer's
retained P1-to-cell increment. Selected-cell release independently checks the
actual native rigid partition, then composes that increment with the complete
selected-piece release. Successive releases can consume retained increments
without subtracting rounded totals. Original inertia, mass/momentum, caloric
state and independent remap budgets remain unchanged; numerical loss never
becomes heat or release work. Four executed production-leaf tests cover hidden
spin/orbital storage, actual mass/momentum defects and invalid inertia. Native
whole/partial regression assertions additionally cover analytic P1 spin,
thermal attachment and successive releases; full native execution is pending.
The source remains an already-liquid no-storage rigid constraint, not a
stressed solid whose constitutive history can be discarded.

## Native solid-electrode circuits with temperature feedback and finite heat

`SpatialMagneticThermalNetwork::step_terminal_circuit_candidate` closes the
original 3D Whitney A–phi electrodes, native RLC graph and finite caloric cells
in one accepted candidate. Every thermal proposal resolves scalar or oriented
tensor conductivity and induced-linear permeability from the original material
sources. Field-response probes retain the same accepted magnetic potential,
native charge/flux and thermal history; their actual count consumes the shared
thermal trial budget. Endpoint operators are prepared once per proposal and
the native field response is invalidated when material coefficients change.

Original local conductor-frame Joule heat enters whole-volume field receivers.
Every authored graph resistor additionally requires one explicit finite heat
destination, which may be a separate thermal body. These physical losses are
deposited once; artificial reference-resistor dissipation is excluded. Native
capacitor/inductor storage and external source work retain their own accounts,
while electrode interface power cancels internally. The existing magnetic
caloric owner preserves signed parameter work, field entropy and reversible
heat instead of treating every magnetic energy change as Joule loss.

The joint receipt checks `delta U_field + delta W_circuit + delta H_caloric
+ W_mechanical = W_impressed + W_external + Q_external`. Original sources,
clocks, volume correspondences, thermal laws, material and solve budgets remain
admitted before joint publication. Refused heat routing, source changes,
cancellation and failed field or enthalpy trials preserve all three histories.

The legacy fixed-resistance entry retains authored R/L/C coefficients,
geometry and conducting occupancy; prescribed motion represents external mechanical work. This is a linear
magnetic, fully solid MQS coupling with endpoint material/enthalpy feedback,
generally first order in the coupled temperature update. Independent Whitney
element elimination tests exercise both electrical property bases, signed
permeability slopes, separate resistor heating and retained capacitor charge.

## Complete heat transport and material resistor feedback for solid electrodes

`spatial_induction::circuit::heat_assembly_candidate` continues the same native
linear A–phi electrode circuit through any original `HeatAssembly`. Every
proposed caloric state resolves contact resistance, integrated conduction
layers, gray-diffuse enclosure radiation and provider-owned finite thermal
links or coolant together with conductor conductivity and magnetic entropy.
The original backward-Euler residual advances every finite receiver once;
there is no independent cooling step or extra temperature history. The full
associated source sample and provider error type remain in the accepted
candidate or typed refusal.

`SpatialTerminalCircuitThermalNetwork::step_linear_material_candidate`
also binds every authored graph resistor to its original electrical material,
basis and whole-volume finite body. The complete heat adapter consumes those
same bindings from its prepared electrical/material network. Accepted graph
resistance and retained source state must agree before stepping. Every trial
pins original claim selection, full fixed conditions, geometry and receiver,
then rebinds the original R(T) without resetting native charge, flux, time or
magnetic potential. A changed resistance invalidates the native response cache.
Source-backed graph resistors occupy distinct caloric bodies outside the field
receivers; their length times area must equal the corresponding original volume.
The earlier fixed-resistance entry retains its original API and numerical
resistances; it refuses to discard an already retained material binding.

Both entry points share the same linear field/circuit driver, graph-source
admission and authored-loss deposit owner used by the original coupled paths.
Conductor-frame Joule loss and actual graph resistor loss enter their explicit
finite receivers once. Artificial interface-reference dissipation is canceled
inside the native circuit and never becomes physical heat. Magnetic parameter
work, entropy and reversible heat retain their original separate accounts.
Finite-link exchange and closed-enclosure radiation are internal; both the
native and physical energy checks exclude the enclosure's numerical net from
external heat. The joint account remains
`delta U_field + delta W_circuit + delta H_caloric + W_mechanical
= W_impressed + W_external + Q_external`.

Fixed geometry, fully solid field receivers, linear magnetic source domains,
reactive coefficients and first-order thermal feedback remain the admitted
scope. Prescribed motion exchanges external work without advancing mechanics.
Thermal boundaries requiring an unowned E/J/B/H coordinate and electrical
interface contacts without their native work owner refuse. Cancellation,
changed source identity or geometry, invalid routing and exhausted budgets
return no partially advanced field, circuit or thermal history.

Focused regressions use independent Whitney-element elimination and two-body
caloric balances for both electrical source bases, continued drive reversal,
R(T), sigma(T), mu(T) and magnetic entropy. Complete-assembly cases exercise
source-controlled contact and emissivity changes, integrated layers, accepted
provider payloads, refreshed finite links, typed source refusal, cancellation
and unowned field-conditioned radiation.

## Physical nonlinear and hysteretic electrodes with finite heat

## Finite-strain Maxwell membrane relaxation into original enthalpy cells

`membrane_thermal::visco::advance_candidate` composes the native
`FiniteViscoMembraneMotion` with the existing finite enthalpy timeline. Each
actual transport, observer or event interval advances the original momentum
equation and Maxwell histories once. Physical per-triangle interval loss is
handed to the matching original caloric receiver before cumulative-ledger
rounding; variable heat capacity, original graded reference mass and volume,
source identity, fully solid phase and the declared temperature interval retain
the existing membrane thermal admission.

The Kelvin--Voigt and Maxwell entry points share the same private receiver,
timeline and heat-handoff implementation; each mechanical owner retains its
native physics. Backward-Euler kinetic loss is reported separately from heat.
Material endpoint defects, pressure path/endpoint work differences and inertial
or equilibrium residuals contribute to an absolute numerical allowance without
cancelling one another. Held supports must have their accepted displacement and
zero velocity throughout the coupled call. Dead loads and admitted current-area
pressure remain held.

The result publishes motion, thermal inventory and observer together. A late
observer refusal, source/inventory mismatch, work refusal or cancellation
preserves all accepted inputs. Mechanical coefficients remain frozen over the
explicitly admitted fully solid interval: this is one-way self-heating, not
temperature-driven coefficient or natural-reference evolution.
`tests/visco_membrane_thermal.rs` covers regional Maxwell heat into variable-Cp
enthalpy, exclusion of numerical kinetic loss, event/continuation history and
whole-call refusal behavior.

`SpatialMagneticThermalNetwork` supplies `step_nonlinear_voltage_candidate` and
`step_nonlinear_current_candidate`. With `magnetic-hysteresis`, the corresponding
`step_hysteretic_voltage_candidate` and `step_hysteretic_current_candidate`
advance actual solid-conductor electrodes, the original 3D A–phi field and
finite caloric inventory together. The shared `SpatialTerminalThermalCandidate`
retains the typed native electrode owner, finite thermal state and full coupled
receipt. Candidate failure or cancellation publishes neither history.

Private read-only history traits let the existing finite-enthalpy drivers read
accepted electrode geometry and sources without exposing isolated field
stepping. Every thermal trial starts from the same accepted magnetic history,
retains electrode charge topology, and resolves scalar/tensor conductivity at
the actual endpoint temperature. Nonlinear electrodes preserve the existing
optional source-separable magnetic entropy law and effective-capacity checks.
Hysteretic electrodes retain vector polarization and require fixed fitted
coefficients over the full admitted caloric interval.

Actual terminal V*I work joins impressed supply once. Native local Joule heat,
and irreversible pinning heat where present, enter their original whole-volume
solid receivers once; numerical magnetic dissipation remains outside deposited
heat. The joint work balance closes magnetic/caloric storage and numerical loss
against electrode work, impressed work and external heat. The nonlinear path
also retains the original reversible magnetic entropy/parameter-work account.
Thermal clocks, complete magnetic source identity, receiver volume, field and
finite-enthalpy budgets all retain their original admission.

The `nonlinear::terminals` and `hysteretic::terminals` modules additionally expose
`voltage_heat_assembly_candidate` and `current_heat_assembly_candidate` for the
complete original contact, layer, radiation and finite-coolant heat provider.
Provider source payloads, typed errors and closed-radiation accounting are
preserved. Electrical graph/contact work without an electrical owner refuses.
These thermal adapters describe stationary, fixed-geometry, fully solid bodies;
material transport, phase transitions and transported magnetic memory remain
outside this scope.

## Native 3D hysteretic winding circuits and finite enthalpy

With `magnetic-hysteresis`, the existing `SpatialWindingThermalNetwork` supplies
`step_hysteretic_candidate` and `step_hysteretic_material_candidate`. The native
3D graph/winding, A–phi/vector-memory field and finite caloric inventory advance
in one candidate. The material form additionally resolves every ordinary graph
resistor from its separate finite body. The complete
`spatial_induction::hysteretic::winding::heat_assembly_candidate` retains original
contact, layer, radiation and coolant source payloads and closed-heat accounting.

The original field thermal residual and the existing winding source/heat routes
serve both excitation forms. Each thermal endpoint resolves actual strand
resistivity, core scalar/tensor conductivity and sourced graph R(T). Magnetic fit
coefficients remain constant over the complete admitted caloric interval.
Neither thermal nor terminal trials reset or partially publish accepted circuit,
field, polarization or caloric history. The shared thermal budget counts nested
field trials, including discarded terminal proposals.

Core eddy, all-tetrahedron pinning, strand copper and graph-resistor heat are
deposited once into their original finite receivers. Insulating dissipative
cores require receivers too. Strand fills must cover their whole original
caloric volume; partial composites require a separate mixture owner. Ordinary
material graph bodies must be distinct from those whole magnetic volumes.
Numerical field/end loss remains outside caloric heating. The stable account is
`delta W_graph+field+end + delta H_caloric + D_numerical = W_external + Q_external`.

All source, original-volume, clock, terminal, finite-enthalpy and total-work gates
precede joint publication. Source changes, cancellation and exhausted budgets
leave every accepted history intact. These adapters remain stationary and fully
solid; they do not infer transported/rotated magnetic memory or a
temperature-dependent hysteresis free energy from a fixed fitted law.

## Three-dimensional hysteretic induction and finite heat

With `magnetic-hysteresis` (also selected by `hysteretic-eddy-heating`),
`SpatialMagneticThermalNetwork::step_hysteretic_candidate` advances the retained
three-dimensional Whitney A–phi field, all vector branch polarizations and
the original finite caloric inventory in one immutable candidate.
`spatial_induction::hysteretic_heat_assembly_candidate` supplies the same solve
with the complete original contact, layer, radiation and finite-coolant heat
assembly. No independent field or thermal history is restarted.

Every endpoint trial resolves scalar or oriented tensor conductivity from its
original material source at the actual proposed temperature. Planar and spatial
hysteretic consumers share the same fitted-source and full caloric-interval
admission. Magnetic branch coefficients remain constant over that supported
interval; changing temperature cannot silently replace the fit or erase its
memory. Complete magnetic source receipts retain the actual endpoint query.

Exact native Joule and irreversible pinning heat enter their assigned finite
cells once. Electrically insulating magnetic tetrahedra also need a real
thermal destination whenever they dissipate. Original whole-volume, fully
solid material correspondence, source identity and clocks are admitted before
thermal trials. A fixed region with unreceived pinning heat refuses.
Numerical endpoint dissipation is a separate energy export and never a heat
source. The stable work account closes
`delta W_magnetic + delta E_caloric + D_numerical = W_impressed + Q_external`.
Both physical heat channels remain individually available in the native step;
`deposited_heat_j()` reports their already-deposited total.

Complete heat providers retain their original typed errors and closed-radiation
balance, and refuse graph electrical work without an electrical owner. Native
field and finite-enthalpy budgets, source failures and cancellation publish
neither accepted history. This is stationary, fixed-geometry, first-order MQS
with source-admitted isotropic vector hysteresis and constant magnetic fits.
Transported/rotated magnetic memory, temperature-dependent fit entropy, phase
changes and empirical device qualification require their respective extensions.

## Fatigue on original finite-rotation elastic rod motion

With `rod-dynamics`, `rod_motion::fatigue::RodMotionFatiguePath` owns one native
rod trajectory and the shared fixed-material-point fatigue observer. Both
centered and eccentric reference inertia retain their original native motion
step, support reactions and energy/momentum accounts. World-frame nodal loads
and prepared material-point spatial loads use their existing kick/drift/kick
owners. The fatigue observer receives the accepted pose, so unforced vibration
still contributes actual signed stress cycles.

Admission observes the supplied mechanical state once, including its preload
and momentum, without reconstructing earlier damage or inserting a zero-stress
sample. The rod observer retains original per-layer axial/bending/torsional
stress, material-frame plane directions, typed source conditions, bounded S-N
rainflow and optional explicitly closed Findley windows. Its explicit
`new_with_transverse_shear` constructor includes the original material owner's
Saint-Venant transverse field for homogeneous isotropic circles and tubes.
It retains both free contours, source roles and whole-layer shear bounds.
Other constructors retain the solid observer's omitted-shear allowance.
No constructor treats a section-averaged force as a local stress or changes
the native shear correction to fit the observation.

A native candidate is accepted only after every gauge accepts its tensor.
Mechanical refusal, late gauge failure, sample/storage budget exhaustion and
cancellation leave pose, momenta, clock, accepted-step count and all histories
unchanged. Ordinary observation clones bounded rainflow residue and rolls back
appended Findley prefixes, without copying the retained tensor window each
step. Clones preserve continuation; provisional endpoint damage does not seal
the live path. Fatigue introduces no mechanical force or heat.

This consumer exposes immutable motion and uses the native conservative elastic
step at fixed source conditions. Viscous constitutive stress, source-temperature
refresh, impulses, contact and fracture are not silently omitted from a coupled
evolution; those operations require their own compatible fatigue continuation
and are not exposed by this path. Caller-selected time steps must resolve the
actual stress extrema and any rotating shear trajectory. A finite plane set
and empirical Miner or Findley estimate do not establish component service life.

## Retained kinetic changes through layer heat and cohesive motion

The rigid third-body shear-layer owner evaluates accepted storage through
`fs_mbd::RigidKineticEnergyChange`, reducing both original bodies before final
rounding. Small layer elasticity and irreversible work survive large common
translation. The finite layer's existing thermal feedback continues to consume
its independent constitutive loss and original caloric inventory; storage does
not define heat. The source-bound scalar, thermoelastic and tensor routes share
this same native mechanical owner.

Deformable cohesive heat, finite adhesive layers, prescribed source-temperature
transfer, opening/cooling, accepted intervals and sustained/event trajectories
consume the native `kinetic_energy_change_j`. They do not reconstruct it from
absolute endpoint totals. Existing bulk/cohesive storage, fracture work, explicit
heat-conversion fractions, source work and caloric loss accounts remain separate.
The original fixed-mass, fixed-reference, source and topology limitations apply.
G3 regression definitions carry actual rigid-layer heat/temperature feedback and
deformable bond-viscosity heat through large unchanged kinetic inventories.
Arithmetic-leaf tests alone do not qualify the complete coupled integrations.

## Native kinetic increments through spatial heat and thermal shells

The spatial viscoelastic thermal owner consumes the accepted native
`ViscoDynamicStep::kinetic_energy_change_j` and
`material_stored_energy_change_j`. These come from the original full or
supported consistent mass and original per-element constitutive increments.
Explicit support and eigenstrain work remain separate, and only independently
integrated Kelvin/Maxwell loss enters the original finite caloric cells.
Unchanged background translation cannot hide a small physical heat transfer
or enlarge the native interval energy allowance.

The native step's read-only `retained_energy_change_j` and `retained_work_j`
carry the combined storage and independently supplied work in double-double
form through the thermal consumer's complete balance. Opposing large kinetic
and material increments are reduced before diagnostic rounding; wrapping their
already-rounded public scalar fields in compensated arithmetic is insufficient.

The reciprocal dynamic laminate owner consumes the original enriched MITC
midpoint's retained kinetic increment directly. Its original director mass,
nodal/director cross terms, material stress work, finite heat cells and launch
coordinate basis remain shared across temperature trials. Reconstructed
absolute kinetic totals cannot reject or erase the original small exchange.
These remain fixed-reference, fixed-mass models with the existing material and
time-resolution obligations. Focused regression definitions cover supported
spatial strain heat and independent director motion beside dominant unchanged
kinetic storage; isolated arithmetic checks do not execute those integrations.

Finite-solid time steps, impulses, force kicks and kinetic drift retain actual
kinetic increments through the unchanged original consistent mass. Endpoint
difference and mean-velocity tails remain compensated through independent
storage and supplied-work checks, including pressure and gravity ports.
Background kinetic storage cannot erase a resolved local exchange. Numerical
integration, impulse-work and force-splitting defects retain their original
gates and never supply physical heat. `FiniteSolidStep::kinetic_change_j`
exposes the same increment used by the original acceptance calculation.

Dynamic source-bound tetrahedral thermoelasticity consumes the native Newmark
owner's retained kinetic increment, independent force/support work and explicit
Rayleigh loss. It retains those quantities through stress-work and combined
caloric acceptance instead of subtracting absolute kinetic diagnostics. Rayleigh
energy stays explicitly exported; local reversible constitutive heat alone
feeds its original material cells. The underlying consistent mass, stationary
support restriction, source ranges and time-refinement obligations are unchanged.
Its finite-gas cavity composition consumes the same retained mechanical storage
and work rather than rounding them again before the combined gas/wall balance.

The original rigid caloric partial-slip owner and its shared-body network now
retain native fixed-inertia kinetic increments across all participating bodies
and through their complete material/thermal balances. They recover relative
midpoint motion from original COM/spin/arm components before common translation
can swallow real slip. Contact free energy, independently computed constitutive
work, per-surface reversible entropy heat and physical microslip heat keep their
original owners; return-map numerical loss never enters a thermal cell. The
existing absolute energy tolerances do not grow with unchanged kinetic storage.
Focused regression definitions compare actual shear/spin impulses, histories,
physical losses and finite caloric endpoints with and without a dominant common
normal translation. Full coupled execution remains separate from arithmetic
leaf checks.

## Original specific heat in coupled capacity admission

Electrical, dielectric, magnetic, piezoelectric, thermoelastic, creep and
tensor/friction-layer thermal admission now use the original caloric interval's
minimum Cp, rather than its average enthalpy/temperature chord. Clipped
temperature envelopes use conservative intersecting-interval source bounds.
Creep retains its outward arithmetic bounds. Existing source identity,
reference mass/density, field-work and mechanical-energy accounts remain
authoritative; no second caloric or apparent field-capacity term is inserted.

String admission with scalar or constant-table expansion, and concentric-layer
thermoelastic admission, evaluate the quadratic TOTAL capacity, including the
original affine Cp slope, at endpoints and its interior stationary point.
Strings with continuous piecewise-affine expansion instead use conservative
bounds over original CTE/caloric interval intersections. These retain the
original expansion slope, native strain envelope and actual modal amplitudes;
immutable bounds are shared across candidates. Actual free strain and entropy
still use the original integrated expansion and current source coefficient.

Local finite-solid and small-strain adiabatic wave gates compare original Cp at
the actual state. Continuous Cp remains differentiable when its temperature
slope changes; an actual capacity jump still refuses a classical wave
linearization.

Native regressions cover continuous source-Cp wave states and a two-knot
capacitor/inductor caloric source whose average capacity hides negative total
capacity elsewhere in its admitted envelope. Refusal leaves the original
thermal and circuit state unchanged. These whole-consumer tests were not
executed in this change because the real asupersync dependency build exceeded
the available memory; numerical-leaf results do not substitute for them.

## Source-dependent finite-solid heat capacity and adiabatic waves

`FiniteSolidMaterial::adiabatic_acoustic_state_with_coefficients` and the
corresponding `ConservativeThermalField` element view retain the original
Neo-Hookean E(T)/nu(T) source derivative owner alongside the existing caloric
chart, density, expansion path and manufacturing reference. The selected
model, complete fixed context, card and actual temperature must match. These
are the same prepared mechanical states used by the ordinary fixed-coefficient
path; the added source capability supplies their thermal derivatives.

The shared native thermoelastic potential combines coefficient curvature,
expansion derivatives and their mixed terms. Its `-T Psi_TT` enters total
fixed-F capacity once, together with the original source capacity at the
actual state and reference density. The full beta=dP/dT and positive total
capacity then feed the existing compressible adiabatic acoustic tensor,
signed modes and energy velocities. Original source receipts remain borrowed
by the result. Nonpositive capacity, foreign evidence, unsupported smooth
neighborhoods and unresolved states refuse.

The original stress-free zero-pressure caloric convention remains explicit.
This is a local reversible full-compressible response; constrained mixed or
incompressible projections, heat transport and attenuation require their
existing separate owners. Finite temperature evolution at held deformation is
available through the caloric transaction below.

## Finite heating and cooling at held deformation

`FiniteSolidMaterial::fixed_deformation_caloric_state` retains one immutable
original material, coefficient source, fully solid caloric curve, total F,
manufacturing reference map, reference mass and reference volume.
`FiniteSolidCaloricState` explicitly uses internal-energy storage with the
existing `CaloricReference::StressFreeZeroPressure` convention:
`U = m*h0(T) + V0*(Psi - T*Psi_T)` and
`dU/dT = m*Cp(T) - V0*T*Psi_TT`. Mechanical caloric coupling enters once.
Each temperature query evaluates the original E(T), nu(T) and natural stretch,
then returns the evolved native stress, tangent and adiabatic response.

`add_heat` accepts signed finite heat and reuses the shared safeguarded caloric
inverse. Original Cp and Cp/T finite-pair integrals cross retained source knots;
direct native mechanical increments avoid cancellation against large absolute
caloric or elastic baselines. A bounded sampled temperature bracket requires
positive combined capacity and increasing energy at every accepted sample.
The final result independently checks heat closure and endpoint Clausius
bounds. Explicit spacing, state-evaluation, root and entropy controls bound
work. Refusal or cancellation leaves the original state unchanged; zero heat
returns that state without claiming an unsampled branch.

This is a fixed-F, fixed-reference-volume transaction on the original fully
solid smooth source branch. Expansion must remain affine across the prepared
path; unsupported kinks, foreign sources, unbracketed heat and nonpositive
sampled capacity refuse. Numerical sampling does not establish global
monotonicity, stability between samples or a thermal CFL bound. Free mechanical
motion, phase change, finite-pressure Cp conversion and spatial heat transport
remain outside this transaction.

Six original consumer tests pass in the native selected finite-solid closure,
using rebuilt original material and thermochemical code with real dependencies.
They cover heating/cooling across a Cp knot, refreshed stress and waves, a
sub-ULP caloric datum, weak mechanical energy increments, constant-source
reduction and atomic source/budget/cancellation refusals. Unrelated FE/time
adapters and thermal-field wrappers were not included in that validation;
this is not a full workspace build.

### Timed reservoir exchange at held deformation

`FiniteSolidCaloricState::exchange_with_reservoir` advances that same original
inventory for an explicit duration against an authored
`fs_thermal::enthalpy::Reservoir`. The descriptor retains its covered area,
prescribed bath temperature and total conductance G in W/K; area does not
multiply G again. One backward-Euler inverse solves
`Delta U(T) + dt*G*(T - T_bath) = 0` with the actual tangent `c_F(T) + dt*G`.
It reuses the finite-heat source cache, bounded branch sampling and direct
energy/entropy increments. Temperature-dependent Cp, elastic coefficients and
expansion continue to affect the same evolved stress and acoustic response.

`FiniteThermoelasticReservoirBalance` retains the actual endpoint power and
`Q = dt*G*(T_bath - T_final)` separately from the independently evaluated
internal-energy change. It exposes the body's original heat balance, endpoint
conduction entropy `Q*(1/T_final - 1/T_bath)`, temporal entropy defect
`Delta S_body - Q/T_final` and combined body-plus-bath entropy change. The
time-discretization contribution is not relabeled as physical generated heat
or boundary dissipation. The existing energy/entropy allowances govern final
acceptance; late cancellation leaves the supplied state unchanged.

The prescribed reservoir has no represented finite heat capacity. Its
temperature may lie outside the material's source interval, while every body
state and the offered branch remain supported. Conductance and bath
temperature are held over the explicit interval; duration refinement remains
independent of nonlinear solver tolerances. Insulation, zero duration and an
equal-temperature bath retain the original state after input admission. This
transaction does not add spatial contact topology or infer a global thermal
slope bound from sampled material capacity.

Three additional original consumer tests cover the constant-Cp backward-Euler
closed form, refinement toward exponential cooling/heating with decreasing
temporal entropy defect, variable-source heating/cooling across a Cp knot,
updated mechanics, exact large-datum invariance, supported body states with an
external bath temperature, no-op behavior and atomic domain/budget/cancellation
refusal. The selected native consumer root passes all nine tests using the
same real dependency and original-source closure described above.

## Circuit-driven spatial hysteresis and finite heat

The `hysteretic-eddy-heating` feature is an alias of `magnetic-hysteresis`;
both select one shared field/heat owner, including complete heat assemblies.
Original fitted magnetic claims are retained at the actual endpoint temperature.
Static applicability and the full caloric source interval are admitted once per
candidate before thermal trials.

With the default-off `magnetic-hysteresis` feature,
`WindingThermalNetwork::step_hysteretic_candidate` and
`step_hysteretic_material_candidate` advance the original native circuit,
vector magnetic history and finite caloric state in one immutable candidate.
The material form also resolves original graph resistors at their receiving
cells' endpoint temperatures. Stranded copper and core conductivity are
resolved independently from their original material cards without replacing
the fitted magnetic sources or resetting remanence.

The existing terminal and enthalpy solvers retain ownership. Every discarded
thermal or terminal trial restarts the same accepted circuit, field and branch
polarizations; nested field work respects the original shared trial budget.
Original cell routes receive core Joule, stranded copper, graph-resistor and
all-triangle pinning heat once. Insulating magnetic cells require finite heat
receivers when pinning dissipates. End and field numerical losses remain
outside caloric deposition.

The combined work check uses the native stable graph/field/end storage
increment, external circuit work and actual finite enthalpy increment.
Source-bound graph resistors require exact original routes, distinct whole
thermal-cell volumes, material selections and physical conductor geometry.
A source refusal, work mismatch, exhausted budget or cancellation returns no
accepted circuit, magnetic or thermal history. Fixed geometry, constant
source-admitted magnetic coefficients and wholly solid caloric charts retain
the same first-order estimated MQS limitations as the field-driven adapter.

## Spatial hysteretic induction on the original finite thermal timeline

The default-off `magnetic-hysteresis` feature adds
`MagneticThermalNetwork::step_hysteretic_eddy_candidate` and
`circuit_enthalpy::eddy::thermal::heat::hysteretic_heat_assembly_candidate`.
They compose the original planar vector-history induction owner with the
existing finite-enthalpy residual and complete heat assembly. Every thermal
trial restarts the same accepted potential and all accepted vector branch
polarizations. Conductivity is resolved from the original material card at the
actual trial temperature; fitted magnetic coefficients retain their original
sources, values and explicit constant-temperature validity. Temperature cannot
silently change the magnetic fit or erase remanence.

Native per-triangle Joule heat and irreversible pinning heat are deposited once
into their actual caloric receivers. Pinning heat from electrically insulating
core cells is included. A dissipating fixed region without a finite thermal
destination refuses the candidate. Original magnetic/thermal cell identities,
whole represented solid volumes, source selections and clocks are checked.
Passive coolant cells can belong to the complete heat assembly without lying
inside the magnetic mesh; existing contact, conduction and radiation owners
continue to resolve at the same endpoint.

The work account uses the native stable magnetic energy increment:
`delta W_magnetic + delta E_thermal + D_numerical = W_impressed + Q_external`.
The report's existing `deposited_joule_j` field includes both physical loss
channels, while `electrical.joule_heat_j` and `electrical.hysteresis_heat_j`
retain their separate values. Backward-Euler numerical dissipation never heats
the material. Complete heat providers retain the original enclosure-balance
gate and refuse independent graph electrical work without its circuit owner.
All source, work-budget, numerical and cancellation refusals publish neither
accepted history. No residual correction or extra empirical core-loss term is
introduced.

This is first-order, fixed-geometry, isotropic estimated MQS coupling with
constant fitted magnetic coefficients. Source-supported electrical transport
can vary with temperature; temperature-dependent hysteresis parameters,
magnetic entropy, phase-changing geometry and measured device validation are
not implied. The focused `hysteretic_eddy_enthalpy` tests use manufactured
scalar-play/P1 balances, physical heat separation and transactional refusals.

## Nonlinear tetrahedral current and finite thermal storage

`circuit_enthalpy::conduction::spatial::response` closes the original
`SpatialConductionResponse` and finite caloric cells at one backward-Euler
endpoint. Each proposed temperature/phase refreshes the same complete pinned
bulk sources, then the native tetrahedral owner closes their actual E/J
coordinates. Boundary voltage/current density and mixed equipotential
voltage/current electrodes retain their native meanings. Immutable geometry,
source identity, query schema, caloric charts and the selected source-coordinate
relaxation policy persist through feedback and continuation.

`step_response_contact_candidate` refreshes both actual contact temperatures,
declared mean temperature and phase coordinates while retaining the accepted
electrical coordinates as native solve seeds. Explicit native nonlinear
contact bindings own any voltage-jump/current-density response. Their uniform
trace admission, source bounds and dimensional constitutive gates remain
mandatory. Each contact requires its original adjacent-volume heat receivers
and an explicit physical partition. Native volume E.J and exact surface work
enter those original finite enthalpies once; constitutive defects are only
diagnostics. Independent terminal source work closes against that deposited
heat before publication.

`step_spatial_conduction_response_candidate` composes any complete
`HeatAssembly`, retaining its source samples, thermal contacts, layers,
radiation, reservoirs and finite links at that same endpoint. Caloric and heat
transport claims must be independent of electrical response/applicability;
this path supplies no reversible field storage. Every native basis solve,
including discarded constitutive and thermal stencils and final fresh
admission, consumes the original common field-work allowance. Source-domain
refusal, work exhaustion, late cancellation and failed balances publish
neither history. Original fixed-geometry and confined-phase restrictions apply.

Planar `ConductionResponseCircuit` retains the complete original bulk/contact
source bundle and relaxation policy across frozen-operator transfers. Its
report retains actual-field contact source evaluations beside bulk source
defects and native circuit work.

## Nonlinear tetrahedral resistance with native circuit storage

`conduction::spatial::circuit::response::SpatialConductionResponseCircuit`
closes the complete original bulk/contact laws at the actual native circuit
midpoint. Every iteration rebuilds the full multi-electrode conductance block,
then proposes the same interval from the same accepted charges and fluxes.
Initial source admission also uses the actual authored electrical state.
Fixed C/L, reciprocal mutual links, original source ports and ideal transformers
retain the existing native descriptor meanings; no pair-resistor replacement
or accepted-resistance freeze is introduced.

`SpatialConductionThermalNetwork::step_response_circuit_candidate` and
`step_response_material_circuit_candidate` place this electrical closure inside
the original finite-enthalpy residual. Complete pinned material sources refresh
at each proposed temperature/phase, and accepted electrical coordinates remain
the source-query seed on continuation. Source-bound graph resistors retain
their own distinct whole-cell receivers. Each surface contact retains both
original temperatures, nonlinear trace ownership and explicit physical heat
partition. Actual volume E.J, exact contact loss and native graph-resistor heat
enter their original finite receivers once. Discarded source/thermal trials
neither advance storage nor deposit heat; source defects are diagnostics.

The linear and nonlinear 3D circuits share one native trial body. Original
terminal-current, factored-power and local-heat gates remain active. Independent
midpoint port work and direct C/L storage increments, including mutual terms,
close the combined electrical/thermal account without subtracting large stored
energy endpoints. Every constitutive iteration reserves N-1 basis solves plus
one circuit step against the common remaining work allowance before starting
that iteration. Source support, finite work limits, nonconvergence and
cancellation refuse the joint candidate atomically.

The bounded lane is scalar, memoryless and fixed geometry. Source-dependent
C/L parameter work, nonlinear tensors, general nonlinear contact quadrature
and electrically dependent caloric/heat-transport laws are not inferred.
Manufactured regressions cover nonlinear RC with two-body thermal feedback,
both E- and J-dependent source laws, polarity and continuation; independent
RL discharge and a three-electrode RLC midpoint; and field/combined work-budget
refusal without modifying accepted storage.

## Relaxed nonlinear electrical source closure

`ConductivitySolve::with_relaxation(weight)` explicitly selects a finite weight
in `(0,1]` for the original circuit and reciprocal field/motion conductivity
iteration. Finite-enthalpy coupling retains that selection through each thermal
trial, its remaining electrical-work budget and final joint admission.
Only bound electrical query coordinates are relaxed. Every next
bulk/interface resistance comes from the same pinned source at that query;
temperature, pressure, source identity, geometry and accepted state stay fixed.
Plain controls and weight one preserve the original arithmetic. Constant laws
still finish in one native trial. Intermediate queries count in
`source_evaluations`, and no query is performed after the native-trial budget
is exhausted. Original actual-field resistance and constitutive-power gates
remain mandatory; native Joule heat is counted once. Source-domain failures,
nonconvergence and cancellation remain transactional. This can resolve passive
laws with oscillating Picard iterates but supplies no general convergence
guarantee. `with_adaptive_relaxation(initial, minimum)` additionally selects
bounded vector Aitken weights through the shared `fs-ports` owner, using all
original bulk and contact electrical coordinate residuals. Only convex source
queries are proposed. Original source support, finite native-trial budgets and
all physical acceptance checks remain authoritative. The finite-enthalpy owner
preserves this policy when narrowing each trial's remaining electrical budget.
Adaptive weights add no electrical solve or source query of their own;
constant laws still finish in one trial. This is an iteration accelerator,
not a stability, uniqueness or time-error certificate. Focused circuit regressions use
strictly monotone terminal laws with independent quadratic midpoint solutions;
the coupled regression deposits their original Joule heat once in a finite
body with source-dependent heat capacity and checks the original enthalpy.

## Finite-solid material potential work

The finite-solid material binding forwards direct elastic energy increments
through the same held natural reference as its energy, stress and tangent,
including the original reference-volume scale. Midpoint candidates use the
prepared solid's single-pass material interval for endpoint storage and its
separately accumulated change. Force-kick/contact drift, mixed free flight and
sustained-contact accounting carry that direct `strain_change_j` into their
original balance equations instead of reconstructing it from rounded totals.
Mixed free-flight kinetic/source/defect sums remain compensated through their
final residual, retaining small work that rounded aggregate totals would lose.
Each original material and total deformation still admits both endpoints.
The change introduces no new constitutive potential, damping, heat or
time-integration scheme; temporal and nonlinear-solve errors remain measured
against the existing physical tolerances.

## Stationary gravitational potential work

`solid_gravity::bind_potential_increment` transports the original native point
field's optional analytic potential difference into the existing tetrahedral
reference-mass quadrature. Uniform and central sources retain small work that
independently rounded absolute potentials would erase. Custom fields without
this capability retain sampled endpoint differences; no additional precision
or gauge-invariance guarantee is inferred for that fallback.

Finite-solid midpoint accounting and force-kick/contact drift now obtain
gravity work from this shared increment. Mixed rigid/solid free flight consumes
the native rigid receipt's potential increment. Endpoint absolute potentials
remain diagnostics. Gravity is still counted once as conservative source work;
the measured difference from the integrator's force work remains a numerical
defect, never heat. The midpoint integrator is unchanged and nonlinear fields
still require temporal and spatial refinement. Endpoint diagnostics no longer
assemble unused Newton tangents.

Reciprocal finite-solid trajectories retain explicit native source snapshots
at both configurations and integrate their direct interaction-potential
increments over the original cellwise or nodal P1 reference-mass quadrature.
Both source recoil and quadrupole orientation participate before rounded
absolute binding-energy totals are formed. Central snapshots retain G M;
quadrupole snapshots retain the original mass, inertia, G, extent and frame.
The native force, recoil, torque, kick/drift stages and domain/energy gates
remain unchanged. The new reciprocal regression covers both mass inventories
and source selections; complete coupled execution still requires the original
constellation dependencies.

The new public-adapter and native free-flight regressions compare zero and
1e20 m uniform datums, including cellwise and graded mass inventories, one-joule
lifting work and 0.125-joule kinetic/potential exchange. Whole-workspace and
full solid/coupling execution remain unavailable with missing constellation
dependencies; these integration regressions are retained for that lane.

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



## Moving rigid damper networks on the finite thermal timeline

`rigid_damping_heat::advance_candidate` composes the native full-spatial
`step_viscous_network_with_loads` with `DepositionSimulation::advance_observed`.
The caller supplies the original body states, complete inertias, material arms,
passive six-coordinate connector laws and held external forces/couples/gravity.
Their accepted time must match the original thermal clock. Every actual thermal
interval, after CFL restriction and exact event/source-transition clipping,
receives exactly one native mechanical step of that same duration.

One authored `BTreeMap<CellId, fraction>` belongs to each original connector.
Fractions are finite, nonnegative and total at most one; compensated admission
retains sub-ulp overspend. Different edges may share finite receivers. Their
independently integrated source heat accumulates once per receiving cell and
becomes average watts over the actual accepted interval. The remainder stays
explicitly unallocated. Neither a physical partition nor a connector-to-material
correspondence is inferred. Existing nonlinear H, fixed-volume U and declared
sensible-storage laws retain their original masses, coordinate bases and caloric
inversions. No instantaneous heat-capacity temperature update or new store is used.

The original observer's source, conductance, ambient, reservoir, pair-transfer,
interval-source, transition, event and timestep hooks remain active. Every mapped
receiver is checked on the actual interval, after all same-time births. Native
authored boundaries and links are restored by their original owner. Each pending
mechanical receipt becomes accepted only in the matching thermal `step` callback.
The native whole-call candidate discards bodies, heat and observer history together
after any mechanical, caloric, observer, conservation, resource or cancellation
refusal. Observer Clone must deeply isolate mutable physical history; external
callback side effects retain the original observer contract.

The raw thermal audit includes generated heat once as subsystem input. The combined
report separates actual edge heat, represented deposited heat, unrouted outflow,
outside thermal supply, original mechanical work, gravity potential storage and
native numerical defects. Per-cell rate conversion and merging with an existing
interval source retain their absolute handoff error. Small physical heat cannot
silently disappear beside a large outside source. Complete mechanical/thermal
energy residuals, accumulated absolute interval mechanical defects and absolute
handoff defects each obey the explicit whole-call joule allowance; numerical
defects never become heat.

The baseline entry keeps mechanical laws, inertias and held loads frozen over
this call. The caller owns their validity interval and subsequent refresh. This is finite
self-heating, not a temperature-feedback damping law, thermomechanical expansion,
phase-dependent mechanics or moving thermal geometry. The existing native thermal
order, geometry splitting and explicit source subdivision still require refinement.

Additional routing work is O(bindings * log(receivers)) per native interval, with
bounded edge arrays and receiver maps. Initial bodies and caller observer remain
borrowed until native candidate work begins; observer physical history is cloned
once, and accepted body banks move from the native receipt. The original thermal
clone and native network allocations remain; no full trajectory, dense global
damping matrix or second thermal solver is constructed. Eight native integration
regressions cover independent exponential losses, shared nonlinear H/U receivers,
partial routing, all observer hooks, event-clipped intervals, native zero-law
motion, same-time births, empty thermal inventory, precision/refusal boundaries
and late-failure/cancellation retry.

## Temperature feedback in moving rigid damper networks

`rigid_damping_heat::temperature::advance_candidate` feeds each explicitly
authored source cell's original caloric temperature into a deterministic,
memoryless `TemperatureDampingLaw`. Each original connector has one declared
CellId; this coordinate is independent of its heat fractions. Shared source
cells are allowed and sampled once per synchronization. There is no average of
two body temperatures, inferred interface temperature or new heat-capacity law.
Original H/U storage, transport, birth/boundary events, observer hooks and exact
native interval selection remain with the same thermal owner.

Each actual interval resolves the complete passive six-coordinate damping form
from the ORIGINAL reference matrix and its accepted start temperature, then
uses the same native moving network and compensated heat handoff. At the complete
endpoint, after all same-time native events, every source resolves successfully
against the actual final caloric temperature before any candidate is returned.
The output includes these endpoint forms and temperatures for inspection. The
caller retains the original reference graph across continuation calls; an
already shifted output is never a new reference calibration.

`temperature::shift::CalibratedKelvinDamping` binds that original graph to one
declared native WLF/Arrhenius `TemperatureShift` per connector. Its original
convention is `a_T = tau(T)/tau(T_ref)`, so the Kelvin resistance is
`C(T) = a_T C_ref`. The adapter checks exact original reference binding, retains
full anisotropy/cross-coupling/null directions, and uses the matrix owner's
retained-factor scaling without another semidefinite factorization. It refuses
source-domain violations, unrepresentable coefficients and compounded shifted
references. Source calibration and the applicability of temperature
superposition remain caller declarations, not inferred material evidence.

This is a FIRST-ORDER lagged-temperature coupling, not an implicit endpoint
solve. A final law consistent with final temperature does not mean that law was
used throughout the preceding interval. Source-freezing error needs independent
max_dt_s refinement even when native mechanical damping and energy conservation
are exact. Native geometry splitting and thermal transport retain their own
orders and controls. Inertias, attachment geometry and applied loads stay fixed;
elastic moduli/reference changes, Maxwell history, expansion, phase-dependent
mechanics and moving thermal geometry require their existing physical owners.

Changing only a memoryless passive damping coefficient stores no energy and
therefore creates no additional parameter work or heat. Original independent
endpoint mechanical/thermal storage, physical edge losses, per-cell handoff
precision and whole-call absolute defect gates remain unchanged. A late source
domain failure, thermal/observer refusal, cancellation or energy failure drops
the entire proposed body/thermal/observer state. Pure source evaluation cannot
mutate external physical history; external callback effects are not reversible.

The source keeps one cached graph and its current per-connector temperatures.
Unchanged temperatures reuse already admitted forms. Additional live storage is
O(connectors + distinct source cells); source sampling is O(connectors *
log(distinct source cells)) per native interval. There are at most
connectors * (accepted intervals + 1) source evaluations and no added mechanical
relaxations, thermal solver, trajectory storage or outer nonlinear iterations.

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
friction-interface heat into both bodies and a solved common interface temperature.
No prescribed friction-heat fraction or arithmetic temperature mean is inserted. The source callback
owns original conductivity claims and other transport, retaining its query and
selection evidence in the returned payload. The route excludes its old direct
link to avoid duplicate conduction. Only actual native Coulomb dissipation
enters the interface heat source. Signed temperature/recession parameter work
remains separate; its explicit external-work budget is the sum
of absolute parameter contributions. Native contact identities and the thermal
cell, interface and aggregate balances have independent acceptance tolerances.

`with_plastic_work_heat(PlasticWorkHeatPartition)` optionally converts declared
fractions `f_A, f_B >= 0`, `f_A + f_B <= 1` of the two native normal plastic-work
increments into heat in the original receiving bodies. The same coupled trial
uses `Q_A = f_A*delta_Wp` and `Q_B = f_B*delta_Wp`; its endpoint temperatures
therefore feed back into original hardness and pressure. These are direct bulk
sources, never extra friction-interface heat. Accepted cumulative plastic work
is not reheated. The callback must exclude both exported friction and configured
plastic heat, while other authored heat at the same cell is retained.
The default remains zero conversion. Fractions are explicit reduced-model
assumptions, not a universal Taylor--Quinney constant, a calibrated material law,
or an inferred split of relative indentation between the two solids. See
Knysh and Korkolis, *Mechanics of Materials* 86 (2015), 71--80,
<https://doi.org/10.1016/j.mechmat.2015.03.006>.

The report separates total incremental normal plastic work, per-body plastic
heat, its nonthermal remainder, and other external heat. The remainder is not
claimed to be recoverable storage. The lower-level thermal report counts direct
plastic sources as external to its thermal-only boundary; `other_heat_j`
removes those sources from the coupled report's external heat. The combined
local identity is `delta_Ucontact + delta_Uthermal + Wp_nonthermal = Wnormal
+ Wparameter_temperature + Wparameter_recession + Wfriction + Qother`.
This retains the externally supplied material-parameter work of a thermal
refresh; its conversion does not establish reciprocal thermoelastic closure.
Focused G0/G1/G3 regressions cover conservative partition admission, asymmetric
body heating, additive authored heat, independent two-cell balances and no
reheating of accepted indentation history.

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
validity. The default heat step does not export material; the optional endpoint
operation below removes actual mass and carried enthalpy. Debris evolution
remains outside this model.
Reciprocal thermoelastic free energy, a resolved subsurface plastic heat field,
moving contact coverage and storage in caller-supplied tangential body response
are not closed by this owner. Its local contact identities and fixed-inventory
heat balance imply no complete body-energy, entropy or passivity theorem.

The candidate seals successor source, contact history and thermal state for
joint adoption. Bounded trials, Newton iterations and thermal inventory,
source-domain or footprint refusal, excessive cumulative removal or parameter
work, nonconvergence and cancellation publish no partial physical history.
The report retains both native mechanical operations, source callback evidence,
friction heat, normal/plastic/parameter accounts and cumulative removed volumes.

### Depth-dependent hardness in live caloric feedback

`friction_heat::plastic_wear::depth::DepthSlidingHeatContact` binds the
original `SourceDepthSlidingWear` history to the same finite-cell caloric
engine. Each trial resolves each body's original H(T,z) at its own temperature
and exposed material depth, with independent elastic properties. Original wear
and friction sources use the solved interface temperature and actual pressure.
Actual unsigned Coulomb slip drives each body's own depth-hardness integral;
signed shear travel, retained indentation and removed inventories survive.

Held approach, maintained force and a fixed finite-stiffness drive use the
native boundary continuations. Temperature refresh, recession and exposure of
the new hardness cap retain three separate normal-work receipts. All three
enter the work audit, sum-of-absolute-parameter-work limit and optional explicit
plastic-heat conversion; parameter work remains separate from friction heat.
The complete thermal, depth, shear and contact candidates are adopted together.
The ordinary and depth owners share the existing Newton/caloric engine and
source/footprint/inventory admission, with no second contact or thermal law.

Fixed and evolving active footprints retain their original geometry limits.
Contact resistance is either explicitly fixed or resolved from the original
interface source as described below. Adaptive native stepping is available for
prescribed actual sliding velocities. This remains a uniform-face model without
lateral temperature fields or a reciprocal thermoelastic free-energy model.
Thermal mass and geometry remain fixed within each interval, and cumulative
removal uses the original receiving volumes. The regression definitions cover
an independent coupled quadratic,
actual reversal under force and spring control, all three work accounts and
transactional source, removal, iteration and cancellation refusal.

### Original thermal resistance on depth-hardness contact

`DepthSlidingHeatContact::with_sourced_contact_resistance` binds the same
canonical resistance source as ordinary plastic sliding, with its original
ordered interface, typed pressure axis, temperature role and fixed conditions.
Each thermal trial queries the source at the final pressure after both native
recession and exposed-hardness continuation. Mean pressure uses the actual
positive-pressure cells on the equal-area grid. The source's original claim and
support remain authoritative; unsupported final pressure or temperature refuses
the complete contact/thermal candidate.

Even with a fixed active footprint, the native heat engine updates both
resistance legs and the successor series conductance. The explicit resistance
split locates the massless generation plane; finite-body temperatures and
conductances determine the heat partition. Accepted steps and material export
retain the same source and split. Ordinary and depth owners share source
admission, pressure sampling and route admission; no second heat law is added.
Four regression definitions cover an independent asymmetric heat quadratic,
actual R(T,p) sampling through wear/export/reversal, zero-resistance replay and
transactional source/support refusal. This remains a uniform-face model without
resolved lateral constriction or third-body caloric storage.

### Adaptive depth-hardness and finite-heat continuation

`DepthSlidingHeatContact::march_prescribed_velocity_candidate` uses the same
native full-step/two-half-step controller as ordinary plastic sliding. The
velocity specifies actual local sliding with zero tangential compliance; the
controller does not subdivide an ambiguous released-body displacement. Only the
two complete native half-step candidates are accepted, preserving original
source identities, both bodies' depth-hardness histories and optional sourced
resistance. Coarse and rejected trials contribute no accepted work or wear.

`DepthSlidingMarchControls` adds independent receiver-temperature and exposed
hardness tolerances to the existing energy, pressure, length, interface
temperature and work controls. The comparison checks each body's original
H(T,z) at its actual recession as well as the effective pressure cap. Relative
scales use changes from the accepted origin, so a large background temperature,
hardness or history cannot hide a small unresolved increment. Temperature
refresh, recession and exposed-cap normal work all enter the accepted span's
contact, thermal and combined energy checks.

Attempt, accepted-substep and minimum-step limits bound refinement. Native
source, nonlinear, resistance and material-domain refusals remain authoritative;
they are not hidden as ordinary accuracy rejections. Cancellation or late
refusal returns no partial span. Interior thermal events require an explicit
boundary, endpoint events remain pending and material export stays an explicit
conservative operation between accepted spans.

Four adaptive regression definitions compare accepted histories with native
steps, refine against an independent continuous graded-hardness/finite-heat
reference, exercise the body-temperature and hardness controls independently,
and preserve predecessors on late refusal. A combined regression continues
adaptive sourced-resistance contact through mass/enthalpy export and reversed
actual sliding, retaining source pressure, route conductance, wear history and
the finite capacity of the remaining bodies.

### Conservative export after depth-hardness sliding

`DepthSlidingHeatContact::export_wear_candidate` removes newly worn material
through the same native receiver-export path as ordinary hardness-limited
sliding. Original sourced density sets removed mass, and the thermal inventory
owner carries the material's current signed enthalpy out with that mass.
Remaining specific enthalpy and temperature stay unchanged, including for a
negative enthalpy datum. This material outflow is never deposited again as heat.

The export checkpoint excludes wear already present at binding and prevents a
second export of the same volume. Original cumulative removal limits survive
export. Explicit shortened prismatic geometry and affected contact patches
retain the ordered bodies, cross sections, contact area and tangent registration;
subsequent heat steps use the remaining masses and new conduction distances.
Hardness profiles, actual recession, indentation and shear history remain with
the original contact owner. Geometry, inventory or cancellation refusal returns
no partial state. Debris evolution and resolved remeshing remain outside this
bounded heat-then-export approximation. Three additional regression definitions
cover signed enthalpy, the next finite-capacity heat balance, repeated/prebinding
export, the cumulative envelope and transactional refusal.

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
interpolation/clamping rules and cancellable validation/evaluation. Canonical
typed scalar payloads use `BoundaryCondition::pressure_at` and the shared
scalar evaluator, including explicit table tails and exact static/instantaneous
Pressure semantics. Dimension-only Pa remains admissible; Stress, AcousticPressure
and amplitude contracts refuse. Scalar pressure is invariant under a checked
frame transform. Spatial profiles require
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

The typed scalar pressure extension adds exact native D3Q19/D3Q27 comparisons
against equivalent legacy histories and checks late source refusal followed by
retry without population, target or clock changes. These added tests remain
unexecuted: the build environment became unavailable during dependency setup,
so Cargo typechecking and runtime verification are still required.

The native `scenario_open_boundary::lattice3` serial/pooled step APIs admit
selected pressure chirps against the frozen `MaterialScaling::timestep_s()`
before preparing boundary targets. `FrameWallVelocity3` performs the same
admission on its prepared ancestor path before writing velocity scratch or
entering the moving-wall operator. The source-owned check covers any interval
of that duration containing the caller's endpoint/midpoint sample, using the
local `sample_time +/- duration` frequency envelope and the default 1/16-cycle
allowance. Chirped tilts also bound individual angular travel; unrelated frame
rows do not constrain the selected wall drive. Typed `FrequencyResolution`
refusals retain the physical sample time/timestep and require explicitly finer
material scaling and a rebuilt grid. Pressure targets, populations, accepted
step counts and wall scratch remain intact on this refusal. Point sampling and
explicit rebinding keep their existing caller-owned integration semantics.
Regressions cover falsely constant pressure, falsely stationary walls, serial
and pooled atomic refusal, and a resolved rebuilt clock with analytic velocity.

The feature-gated `scenario_pressure::compressible` adapter also resolves
native chirp pressure and temperature histories before requesting an Euler
reservoir interval. Both selected ends constrain the shared physical clock to
the scenario default of at most 1/16 cycle per interval, while existing sweep
joins and exact held-value left limits remain in force. Native CFL and retry
logic may shorten that interval further; actual RK stage times still drive
the original source sampling and reservoir flux. Source-resolution failures
retain their original signal diagnostic and publish no fluid state. The same
private advance owner serves the nonlinear duct adapter. This is authored
source resolution, not a general nonlinear bandwidth or time-error bound.
The native reservoir regression checks independent pressure/temperature
cosines at actual stage times together with mass and energy conservation.

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

## Bandlimited fractional latency across pressure ensembles

`PressureEnsemble::new_with_bandlimited_source_delays` aligns existing numerical
receiver latencies that become fractional on a common output clock. For example,
the impact observer's 80 source-sample derivative/interpolation latency becomes
80/3 or 80/6 output samples when the original source runs at three or six times
that clock. The caller declares a resolved physical pressure band at or below
0.45 of the output rate; this declaration does not remove higher source content
or prove convergence. Original decimators retain their antialiasing histories,
source-group phase, mechanical clocks and event order.

Fractional alignment pads each existing numerical latency to the next integral
sample at or after their maximum, then retains a common 32-sample reconstruction
delay. Physical source-to-receiver travel is excluded from that alignment and
keeps its original phase. The original integer constructors and all-integral
requests keep their original integer-delay path bit for bit, without added FIR
latency. Existing whole-sample metadata is accompanied by explicit source and
alignment fractions and the separate reconstruction delay. Both final numerical
latency and complete per-part history remain within the 8192-sample limit.

Physical impact propagation and ensemble alignment share one prepared
`pcm_wav::delay::FractionalDelay` implementation using the existing Rayleigh
kernel. Preparation is cold; callbacks allocate nothing in the alignment layer.
The compensated multi-source sum still preserves weak signed pressures at
destructive interference. Original finite horizons, unflushed tails and failed
callback poisoning remain in force.

Five new native regressions passed using analytic monopole receiver potentials
through the original derivative, physical propagation, 3x/6x decimators and
ensemble runtime. Relative waveform error stayed below 1.501e-4 through 0.45 of
the output rate with the full 139-sample common numerical delay retained.
They also cover an integral member sharing fractional reconstruction, callback
partition equality, source sample accounting, a 2^-60 pressure between cancelling
sources, integer-path parity, band refusal and complete-history limits. All 23
original shared decimator/path/delay tests passed. This bounded native build
retained the actual numerical/runtime owners and exercised original error
variants, omitting unrelated voice facades and error payloads; it is not a full
Cargo or mechanical integration gate.

## Bandlimited harmonic pressure from impacts

The solver-neutral spherical-harmonic radiation bridge fits each channel
through `fs_vfit::vector_fit_transfer`. Observation coefficients retain either
sign, including direct feedthrough; the impedance-only nonnegative direct-term
constraints do not apply. Stable poles, properness, declared direct-term policy
and the independent complex-response gate remain authoritative. Opposite modal
orientations preserve the field when their source coordinates also reverse,
and equal drives on opposite radiation rows cancel.

`render::plate::impact::audio::ImpactPressureRenderer` applies the original
acceleration-response radiation bank to accepted modal velocities, sums its
fixed spherical-harmonic receiver weights, and differentiates the resulting
potential in Pa s with `pcm_wav::decimate::Decimator::new_differentiating`.
The immutable borrowed bank, fixed projections and stationary receiver make
this commutation valid for the whole run. One derivative serves all source
and radiation channels; the callback retains its prepared buffers.

The shared derivative replaces the backward difference's frequency-dependent
sinc loss. Its 48-mechanical-sample group delay is reported by the existing
`velocity_difference_delay_s` timing field, separately from source-radius shift,
geometrical travel and propagation reconstruction. Fitted bands beyond 0.45 of the
mechanical sample rate refuse. Source harmonics already aliased by mechanics
remain unresolved; source-clock refinement and anti-aliasing before a rate
change are still required. Mechanical/contact/thermal histories and energy
accounts stay with their original owners. This frozen-reference acoustic model
does not update its transfer from evolving geometry or material properties;
time-varying transfers require their own derivative product terms and are not
admitted by this interface.

The receiver retains the original `(range - source_radius) / sound_speed`
travel time rather than rounding arrival up to a whole mechanical sample.
Its fixed delay line reuses the existing 65-tap Rayleigh path interpolation
kernel, whose complex-response error is below 4e-4 through 0.45 of the clock.
The additional 32-sample numerical delay is retained even for integral travel,
so different receivers share a common reconstruction latency. Timing reports
whole and fractional physical samples and `propagation_interpolation_delay_s`
separately; `propagation_rounding_s` is zero. The historical nonnegative
`maximum_delay_error_s` bound therefore accepts zero at fractional distances.
`maximum_delay_samples` admits the complete travel/FIR history before allocation.
One prepared mono ring serves the fixed summed field, with no per-mode delay
filters, per-sample allocation or future samples. The physical pressure ceiling
is checked after propagation as well as after differentiation.

Focused regressions use the analytic compact-monopole Green function and the
original exact linear-impact source to check amplitude and phase through the
upper admitted harmonic band, unchanged mechanical states and energy records,
vanishing pressure under coasting motion, fractional listener distances,
separate numerical latency and explicit band refusal. Native propagation tests
also retain the signed interference between paths whose fractional travel
differs even when their rounded integer delays would coincide.

Bounded native verification executed four actual radiation-bank regressions
and two real-bank/shared-derivative component tests. Relative pressure waveform
error remained below 1.54e-6 through 0.45 of the sample rate with the declared
delay; settled coasting pressure was zero. The full impact/contact/thermal
integration tests were added but could not run within the available shared
storage. This is scoped acoustic-owner verification, not a full Cargo gate.

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


## Native nodal heat and continuing tetrahedral dynamics

`solid_enthalpy::spatial::SpatialSolidEnthalpy::resolve_at` binds an explicitly
timed native nodal enthalpy field to the existing source-backed solid state.
`SpatialSolidState::advance_transient_candidate` passes that state to the
original tetrahedral Newmark thermal-update owner. Displacement, velocity,
acceleration and step history continue; the original operator owner retains
reference mass and rebuilds changed stiffness.

`SpatialSolidState::reduce_pressure_valve` also carries the native heat owner's
retained P1 reference-density profile into the common thermal pressure reducer.
Its exact consistent mass controls the pressure response and existing
`SpatialHeatedValveMotion` state transfer; original nodal caloric mass controls
heating. A constitutive density query at the new temperature changes neither
inventory. The G1 graded-inventory definition in
`tests/spatial_solid_enthalpy/motion.rs` compares both owners against an
independent heated one-DOF oscillator and checks late refusal atomically.

The previous mechanical material/eigenstrain state and clock must match the
accepted run. Both endpoints retain their reference mesh, frame, material
assignment, caloric chart, source laws, strain envelope and reference inventory.
Untimed equilibrium states, stale history and incompatible replacement states
refuse. Native support, load, cancellation and energy checks remain active;
failed proposals publish no heat or motion state.

This is prescribed one-way thermal loading on stationary reference geometry.
Signed thermal/stiffness work and Rayleigh loss retain the native mechanical
accounts; none is silently deposited into or removed from nodal enthalpy.
The caller accepts heat and motion together and resolves thermal events and
time refinement explicitly. It does not add reciprocal thermoelastic feedback
or a moving-mesh heat law.

`SpatialSolidEnthalpy::advance_transient` completes that handoff as one candidate:
it admits the caller's accepted nodal heat against the retained mechanical fields,
solves the original heterogeneous conduction/radiation owner, resolves its actual
thermal endpoint and continues native motion using exactly the same supplied
timestep. It never reconstructs that duration by subtracting rounded clocks.
All source-bound heat consumers share the original conductivity/inventory binding.

`SpatialSolidTransientCandidate` returns the original `EnthalpyHeatStep`, timed
source state, native `TetTransientCheckpoint3` and original mechanical work report
together. The checkpoint retains the trajectory and existing shared matrices
across caller-owned material lifetimes; rebinding its exact endpoint problem
continues the accepted clock and step budget. A late structural or cancellation
refusal publishes neither heat nor motion. Original radiative energy reports and
mechanical work remain separate. Regressions cover exact-duration clock rounding,
two-step continuation/replay and structural refusal after the heat solve.

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

For the shared native D3Q19/D3Q27 owner,
`scenario_open_boundary::lattice3::{step_serial, step_pooled}` samples the same
typed histories and advances the original fluid operator as one transaction.
Candidate targets drive collision/stream reconstruction and measured reservoir
transport. Held controls, populations and completed-step count publish together
only after the final native checkpoint; a late local-Mach refusal or cancellation
retains the entire accepted prefix. `BoundaryStep3` returns the accepted physical
samples alongside the original native mass/momentum transport receipt.

All generic mixed and pressure-only adapters reject a selected installed pressure
profile, even when a sampled uniform target equals its dormant default. Unselected
profiles keep their original face fields and still drive the step. The caller
retains the physical sample time and material scaling; no pressure-work or thermal
energy certificate is inferred. Focused tests exercise changing histories on both
lattices, serial/pooled parity, late refusal, cancellation/retry and profile
ownership across the mixed and pressure-only entry points.



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

## Pressure-dependent resistance in hardness-coupled friction heat

`FrictionHeatRoute::from_link_with_resistance` admits two finite nonnegative
area-specific resistances around the existing massless heat-generation plane.
Each thermal leg is `A / (d_i/k_i + R_i)`; the accepted direct link retains the
same series sum `R_A + R_B`. The original direct path is suppressed once during
the native heat solve. Zero resistance preserves the original half-cell
arithmetic. The declared resistance split locates the generation plane; it does
not prescribe the resulting heat fraction or add an interface heat store.

`PlasticSlidingHeatContact::with_sourced_contact_resistance` pins the canonical
thermal-resistance claim from the original ordered interface, original fixed
conditions and explicitly typed static-pressure axis. Every thermal trial uses
the solved interface temperature and mean pressure over the actual positive-
pressure cells of the worn native endpoint. Equal native cell areas make this
mean exactly the resultant/loaded-area convention, without a nominal-area or
pressure/hardness area substitution. Sources calibrated to different pressure
or temperature conventions are outside this reduced model.

The same native cell-energy and massless-interface residuals own all transfers.
Accepted reports retain the actual resistance query and split; the successor
thermal patch and route retain its series conductance. Adaptive continuation and
actual material export preserve this resistance/source checkpoint. Claim support,
identity, geometry or nonlinear refusal publishes no partial physical state.
The model still has two uniform bulk receivers, beginning-pressure sliding/wear,
endpoint heat coefficients and no resolved lateral constriction or third-body
storage. Physical time/grid refinement remains necessary. Focused regressions
cover an independent asymmetric quadratic heat oracle, actual worn pressure,
export/continuation, exact zero-resistance replay and late source refusal.

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

The spatial midpoint owner evaluates kinetic increments through its unchanged
full consistent mass and the shared `fs-math::quadratic` kernel, including all
prescribed-support and free-coordinate cross terms. Elastic/electric increments
use the original local `C^E` and `epsilon^S` tensors at the recovered endpoint
fields. Native Maxwell/Debye storage increments survive integration into the
cell account; their independently evaluated physical heat is unchanged.
Applied-force, support-reaction and electrode work retain endpoint low words.
The ordinary, resistive, reactive and finite-caloric gates consume those
independent increments before scalar diagnostic rounding. They neither subtract
large absolute kinetic inventories nor widen tolerances with background motion.
External capacitor/winding storage remains separate from conductor heat, and
the shared internal thermal-driver boundary retains storage minus source work
until the caloric balance is complete. All original field, charge, source,
amplitude, cancellation and publication checks remain authoritative.

Focused G0 definitions in the private `transient/energy.rs` exercise coupled
tensor cancellation, endpoint work and a small storage/heat exchange. G3
definitions in the existing support-motion suites carry a small transverse
kinetic change beside large common rigid motion through native sensing,
resistor/capacitor harvesting, intrinsic ceramic loss and actual finite caloric
cells. These synthetic dynamic-range fixtures do not establish material
calibration, relativistic physics, integration error bounds or performance.

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

## Prescribed-temperature moisture histories

`MoistureSlab::at_prescribed_temperature` and
`MoistureCartesian::at_prescribed_temperature` propose a uniform-temperature
event at the accepted clock. They retain every cell's dry-basis water content,
dry skeleton, reference geometry, time, definition identity, and original
density/stress-free datum. They re-evaluate the originally pinned diffusivity
claims at the new absolute temperature and every occupied moisture state;
Cartesian fields retain all three explicitly selected principal sources. Any
unsupported cell, source-domain hole at an occupied endpoint, invalid
temperature, or inadmissible diffusivity refuses without publishing a state.
The transport anchor is an occupied cell, so the historical reference moisture
does not unnecessarily constrain a new temperature cut of the transport source.
Reference density is never queried to create or remove either inventory.

`temperature_kelvin` reports the current environment; `transport_context`
exposes its pinned source bundle, while `reference` remains the original datum.
Ordinary and SDIRK stages, finite-atmosphere/solution participants, paired and
adaptive reservoir controllers, layered sorption admission, and sorption-energy
admission all use the current context. Existing boundary laws and reservoirs
must independently qualify the same temperature. Kirchhoff paths are prepared
from the current source context on the next proposal; no old-temperature chart
is retained. Return temperature events use the same original claim pins.

Each transport interval remains isothermal: callers cut intervals at prescribed
temperature changes. This API supplies neither a solid heat equation nor the
heat needed for the imposed event, and does not interpolate an environmental
history across an unqualified source path. Native elasticity can resolve its
qualified current environment; the existing fixed-temperature hygroscopic
expansion integral still refuses a temperature change from its stress-free
datum. Thermal eigenstrain is never inferred from a moisture coefficient.
Authored G1 regressions exercise analytic hot-source flux, principal-axis
refresh, degenerate-front wetting, history invariance, late-cell refusal,
boundary/reservoir temperature agreement and the retained strain refusal.

## Second-order slab diffusion with finite vapor

`FiniteMoistureAtmosphere::advance_sdirk2` couples the original one-dimensional
slab and its actual finite vapor inventory at both Alexander SDIRK2 stages.
The native O(N) tridiagonal diffusion solve, pinned source D(M), optional
Kirchhoff integral, sorption and finite surface-transfer laws remain in use.
Each stage resolves its own RH and material state. Stage two retains the
original solid and gas storage and adds the actual first-stage face load once;
the intermediate gas/solid reports are diagnostic, not additional inventories.

The first-stage physical moisture envelope remains available to the second
stage's source integral. Only current boundaries with nonzero conductance
extend the range check. Invalid stages, source trials, work limits and
cancellation refuse the whole proposal without clipping, partial publication
or a first-order fallback. The returned slab keeps its original ordinary
integration mode. Temperature, gas volume and reference skeleton remain fixed;
sorption heat, condensation and a global temporal-error certificate are absent.

`advance_shared_sdirk2` applies the same native slab stages to multiple
specimens and one common RH, retaining each original geometry, source, signed
face transfer and ordinary integration mode. Common nonlinear controls apply
to every specimen; aggregate cell limits are checked before cloning.
Only first-stage rates and receipts remain during the second scalar search;
diagnostic first-stage fields are released.

`advance_sdirk2_adaptive_to` and `advance_shared_sdirk2_adaptive_to` use the
existing second-order temporal controller. Every cell, RH and vapor mass has
an independent local error scale, so opposed uptake and drying cannot hide
behind an unchanged atmosphere. Only accepted physical half intervals enter
the cumulative signed face accounts. Both implicit stages are charged on
every probe; one fixed whole-duration water, carrier and total-mass budget
admits all endpoints together. A late work limit or refusal discards the
entire proposed duration. The error estimate remains local and asymptotic.
G1/G3/G4 authored fixtures additionally specify heterogeneous shared slabs,
coupled mean/difference modes, stationary-vapor exchange and late refusal in
`tests/moisture/atmosphere/shared.rs` and `tests/moisture/adaptive/slab_sdirk.rs`.

Added regressions in `tests/moisture/atmosphere.rs` specify exact coupled
two-store uptake/drying, second-order refinement, weighted face/vapor accounts,
Kirchhoff bracket support, replay and transactional refusal. These authored
numerical fixtures have not been executed in the editing environment and do
not establish experimental validation.

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

### Physical momentum in reciprocal electric and magnetic gaps

`material_circuit::motion` confines nonlinear field work to the original
charge/linkage and gap coordinates. Its cubic uniform-gap energy uses an exact
average-vector-field chord with an analytic endpoint derivative. The mechanical
effort is the mean physical momentum divided by the original body/support mass;
unrelated circuit storage retains its own quadratic efforts. A full-state
Gonzalez correction no longer changes mechanical velocity to close field energy.
For a held current or voltage source, the actual force impulse matches the
analytic integral of the original quadratic Maxwell force, including reversal.
`motion::dielectric` uses the same average-vector-field rule on each original
piecewise-affine E(D) segment. Exact polynomial quadrature splits only at source
knots; kinetic effort and unrelated nonlinear capacitor secants retain their
physical coordinates. Its linear material limit matches the uniform-gap owner.

`motion::coupled_gap` uses symmetric fixed-gap/fixed-linkage work chords through
the same winding leakage and optional reversible-core source. Core saturation
work is projected only through its original scalar flux constraint. Fixed-linkage
gap chords remain attractive through field polarity reversals. The no-core
case has an analytic endpoint action; saturation and nonlinear dielectric paths
retain the native bounded Newton probes. All four shared-gap corner states must
remain inside original source support; an unsupported corner refuses the trial.
These symmetric, consistent time maps retain the original support, circuit
constraint and transactional refusal boundaries and require time refinement.
Focused regressions check physical displacement/momentum kinematics, reciprocal
energy and source work, polarity and source-knee crossings, and the existing
independent RK4 refinement reference. These are numerical-law checks, not new
measured material or actuator validation claims.

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
Both pre-source and final admission use that envelope. Kirchhoff source charts
also retain the admitted first-stage envelope during second-stage finite-vapor
and finite-solution probes. The carried first-stage flux can legitimately move
a trial outside its current reservoir envelope; it does not authorize source
extrapolation or relax physical stage bounds. Inadmissible stages refuse
without clipping or committing a partial atmosphere/specimen state.
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
1.5e-15 relative error. The five scalar source-to-native regressions and four
shared-core regressions pass in a focused build of their actual source and
numerical owners. This does not establish the complete Cargo/DSR gate.


### Shared-core nonlinear mutual inductance

`material_circuit::inductance::MaterialSharedCoreInductor` binds multiple
original winding branches to one reversible material path. Signed turns and
positive explicit leakage obey lambda_i=N_i*Phi+L_leak_i*I_i, with
sum(N_i*I_i)=l*H(Phi/A)+R_gap*Phi. The existing one-turn source mapper supplies
the exact signed Phi/MMF knots and field-only source admission; the native
`fs-phs::circuit::inductance::SharedCoreInductance` owns their coupled potential.
Saturation from any winding changes every self/mutual differential coefficient.
Opposed ampere-turns can cancel common-core magnetization while retaining both
windings' leakage energy. Source support limits the common flux, not individual
total linkages carrying differential leakage current.

`ReversibleInductorCircuit::with_shared_cores` keeps the original branch
identities, conserved flux states, source ports, Kirchhoff constraints and
native time stepper. Scalar inductor/capacitor curves, reciprocal capacitance
blocks and disjoint linear mutual groups retain their existing owners.
Common-core observations expose B/H, simultaneous currents, the reciprocal
differential inductance matrix, one core/gap energy and per-winding leakage
energies. These partition the native Hamiltonian rather than add stores to it.

All winding ownership must be disjoint. Zero leakage, source-domain exit,
unsupported source dependence, overlapping blocks, external mutual links on
bound windings and nonrepresentable constitutive arithmetic refuse explicitly.
The reduction is fixed-geometry and reversible; it does not infer iron loss,
hysteresis, remanence, fringing, thermal magnetic entropy or a full field model.
Existing moving-gap and hysteretic owners remain separate physical models.
Focused regressions cover source-to-core energy partition, the one-winding
limit, a loaded transformer crossing saturation with native source/heat balance,
and source/ownership/cancellation refusals.


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

## Finite-body thermal feedback through spherical Hertz spin events

`friction_heat::sphere_plane::propose_sphere_plane_hertz_heat` composes the
original sphere Hertz support trajectory with the original finite enthalpy
network. Actual spin arrest, static retention, breakaway, motor reversal and
attitude remain with the native mechanical owner. Every coupled trial starts
from the same opaque accepted motion. Its solved common face temperature
refreshes the original pinned static/kinetic claims once; that snapshot retains
the actual native normal reaction and separately declared footprint radius.

Only integrated sustained physical friction loss enters the existing massless
face. Both finite bodies retain their own temperatures and conductivities;
signed body heat includes ordinary hot-to-cold exchange. Native motor work,
external thermal sources and the independently bounded static constraint-work
residual remain separate. The latter is never thermalized. The combined gate
reconstructs actual compensated native kinetic change plus caloric change,
less outside heat, mechanical supply and the explicit static residual. Motion,
thermal cells, original source receipt and material callback payload leave the
owner together after all original and combined gates pass.

This is a stationary circular contact-column, area-averaged thermal model.
The original thermal face must match the actual Hertz centre, normal and radius;
ordered thermal body IDs and the sphere-A/plane-B side geometry are checked.
The exact pure normal-axis spin subspace admits no COM or plane motion,
tangential drive, impact heat or displaced footprint. No velocity, force or
momentum is projected to meet these limits. Original material/caloric identity
and the applicability of the contact-column approximation remain with their
material owner. Radial flash gradients, material transport through a spinning
footprint, normal compliance, thermal expansion and whole-sphere thermal
geometry are not inferred. Endpoint coefficients are held over each native
window; backward-Euler thermal transport and thermal sampling require temporal
refinement independently of exact held-coefficient stopping.

Five focused G1/G3/G4 native tests pass in
`tests/sphere_plane_friction_heat.rs`: independent analytic spin/attitude and
two-body conduction, first-order convergence against an independently integrated
continuous thermal-fade ODE, heat-induced static breakaway with motor work,
opposed-motor reversal and signed-axis covariance, and late source/cancellation
rollback with exact retry. The 8/16/32-window maximum physical-state errors are
0.0531229, 0.0268324 and 0.0134851. These verify the bounded model with synthetic
source cards, not measured calibration, full-workspace Cargo/DSR or throughput.

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

The optional reed-card record `open_hole_orifice N` (N = 1..=16, `fs_duct::exterior::MAX_ORIFICE_PORTS`; after `constriction_noise`, before `compile_limits`; once
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

An absent row and `N = 1` use the same first-hole shunt, bit for bit (`RegisterHoleShunt::Single`), and an N above
the number of open holes uses all of them. Both the scalar and matrix realizations bound the complete reflected
wave: the direct operator is clamped to power-wave norm 0.999 and the whole tail scaled by the largest factor (found by
bisection, 2^-14 resolution) for which `max_w || D + lambda T(w) ||` stays at most 0.999 on the 2n grid. This is a
sampled realization check, not a continuous-band passivity certificate or a physical material-loss law. The former
first-hole rule, `|D| + max_w |T(w)| <= 0.999`, ignored phase interference and could damp an already passive
reflection. For example `H(z)=0.5+0.5*z^-1-0.3*z^-2` has exact maximum gain `sqrt(68/75)<0.999`, but that rule retained
only 0.62375 of its delayed terms. The phase-aware rule keeps this impulse response unchanged and still limits a
truly active positive-tap response. The scalar specialization shares the matrix path's factor search and stores
only scalar complex spectra; it never allocates a maximum-sized port matrix at each frequency. The old rule remains
`PassivityRule::Triangle` for regression comparisons only. The existing nonlinear single-vs-one-port test now uses
the phase-aware rule on both sides. This changes the single-hole numerical realization while preserving the same
physical TMM, orifice law, explicit loss parameters and scalar runtime. The row takes a constant-medium run or a `gas_field continuous`
binding (`GasFieldDuct::response_with_wall_holes_and_ports`: each orifice reads its own local gas, K and numerical
reference included). A molecular bore, a midpoint `gas_field` (N > 1) and a segment-gas duct with several regions
refuse by name. On a continuous field with `open_hole_interaction` the ports sit on the COUPLED network (`Y` is
`fs_duct::interaction::port_admittance`, the network the flows and the reed's `s12` come from; the independent `Y` is kept only at
a frequency where the coupled network is singular, and at a frequency where the coupled solve is active, i.e. returns a negative input
resistance, where the independent holes are solved instead: bead `frankensim-1xe5t`, listed on the observer by
`OutletPressure::interaction_fallback_hz` and printed in the provenance as `outlet_receiver.open_hole_interaction_fallback`
`{bead, count, hz}`); without the row `Y` is the independent network of the chimneys. Cost: `N^2 + 2N` kernels of the grid length instead of 3
(15 at N = 3, 143 at N = 11, 288 at N = 16), and one more frequency sweep at build; the passivity step holds one `NP x NP` complex matrix
per frequency bin (4.1 KB at the cap of 16, about 67 MB on a 8192-point grid, about 540 MB on the largest 65536 grid; 2.3 KB at the cap of 12 that
preceded it). The cap is a budget (stack arrays and that matrix), not a validity limit. It was 3 until the oboe D5 page (11 open holes), then 12 until the
clarinet G4 (14 open holes) and Gb4 (13); the solver residual test runs n up to the cap
(1200 random passive systems) and `compose::port_tests` checks an eight-hole and a fourteen-hole ladder (reciprocal, passive, leading block
independent of the port count). Measured on the oboe D5 card (11 open holes, `gas_field continuous`, tmmi 2, 96 kHz, 2.65 s):
`open_hole_orifice 3` renders the posted wav byte for byte on the build with the cap at 12, and `open_hole_orifice 11` (all
open holes) changes the steady tone by +0.8 c, the level by 0.00 dB and the harmonic levels by none to the digit, with the
grader count unchanged (6 of 14); user CPU 591 s at N = 3 and 1118 s at N = 11, so the cost is far below proportional to the
`N^2 + 2N` kernel count (15 to 143). Measured on the clarinet pages (24-hole Buffet bore, `gas_field continuous`, tmmi 2, 48 kHz, 2.65 s; 13 pages
with N from 2 to 12 equal to their own open holes, the controls byte-identical to the posted pages): user CPU 4.9 to 5.6 times the control's at every N
(3330 s at N = 3 to 4158 s at N = 12, controls about 700 s; Bb4 at N = 2 took 6609 s against a control of 1340 s), peak RSS 61 to 208 MB; the cost does not
scale with the audio length (A5 and Bb5 rendered for 0.9 s took 3551 to 3693 s at N = 2 to 6, the same as for 2.65 s) or with N, so it is a fixed cost of the render's build, not of the per-sample solve;
the steady tone moved by 0.5 to 4.4 c and the level by 0.2 to 5.5 dB on eleven pages, and on A5 (N = 10) and Bb5 (N = 11) the card played a low tone, 325 Hz 8 dB over A5's own line and 339 Hz with no line left at
Bb5 (the cause is not isolated; the family card records the same kind of tone for F6 without this row). Not measured: an N above 12 on a real card, and an N above 3 on any other family.

Evidence (G0/G1): `fs-duct` `compose::port_tests` (the ladder's admittance is the inverse of the one-port network
impedance to 1e-9, reciprocal, passive, independent of the port count, and the sign of its off-diagonal term fixed by a
coincident-pair test); `orifice_ports::tests` (the solver's residual over 1200 random passive systems, n = 1 up to the cap, and a stiff drive;
the realized matrix is passive and is the admittance scattering, not its transpose; the time-domain steady state with
linear drops equals the network built from the realized kernels to 1e-14 and the TMM chain to a few percent at a 4096
grid; the general solver run as one port follows the first-hole line through a strongly nonlinear drive to 1e-6);
card level `reed::orifice_tests` (an absent row and `open_hole_orifice 1` are bit-identical in audio and reed state, a
one-hole duct ignores the count, 2 and 3 differ from 1 and from each other, bad counts refuse, a fourteen-hole card takes fourteen ports); on a continuous field
`reed::orifice_tests::on_a_field` (an absent row and `open_hole_orifice 1` are bit-identical, 2 and 3 act and differ,
3 composes with `open_hole_interaction tmmi` and the printed lengths, a midpoint `gas_field` refuses by name) and
`fs-duct` `profile::gas_field::ports_tests` (uniform field = homogeneous ports bit for bit, one port = inverse of the
network impedance on a graded field, three ports reciprocal, passive and independent of the port count, the interaction
leaves `Y` alone and changes the radiated flows).

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

Native slabs use `advance_sdirk2` and `advance_shared_sdirk2` through the
same finite-solution root and existing tridiagonal diffusion owner. Each
specimen retains its original material, geometry, surface laws and ordinary
integration mode, including optional Kirchhoff transport. Aggregate cells are
bounded before private staging; full diagnostic first-stage fields are released
after their physical loads and native receipts are extracted. The same
algebraic-phase admission and explicit finite-rate kinetic refusal apply.

`advance_sdirk2_adaptive_to` and `advance_shared_sdirk2_adaptive_to` bind
those native slabs to the existing second-order controller. Every cell,
solvent mass and activity has an independent local error scale, including
pure solvent whose activity stays exactly one. Only accepted physical half
intervals enter the cumulative signed water accounts; both stages are charged
for every attempted proposal. A late work limit or refusal discards the whole
duration. Added authored fixtures in `tests/moisture/solution/slab_sdirk.rs`
and `tests/moisture/adaptive/solution_slab_sdirk.rs` specify pure/phase-buffered
limits, original nonideal activity/osmotic response, second-order refinement,
opposing specimens and transactional refusal; they do not supply experimental
validation or a global/phase-crossing temporal guarantee.

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

## Diffusion between two finite solution baths

`FiniteMoistureSolution::advance_paired`, `advance_paired_cartesian` and
`advance_paired_radial` couple one original native moisture field to two
independently evolving aqueous inventories. Each face names its bath explicitly;
its signed native transfer debits only that bath. A stationary solid-water total
can therefore retain nonzero through-flow between solutions with different
activities. Each bath keeps its own ideal/nonideal activity law, original
nonvolatile stocks and optional equilibrium or finite-rate phases.

The original safeguarded solvent-volume root and phase audit are shared with
the existing single/shared-bath paths. An inner first-bath root is resolved for
each outer second-bath trial; every trial restarts all original inventories.
The two independent solvent balances, field/face balance and compensated direct
whole-system endpoint/start water balance must all pass before either bath or
the specimen publishes. One cumulative field-evaluation budget covers every
nested trial, including rejected ones. Positive source-supported brackets are
never widened; an inadmissible or unbracketed inner probe can refuse even when
a coupled root exists elsewhere. Sealed baths retain solvent while their own
finite-rate phase kinetics may continue.

Native slab/Cartesian/radial diffusion, actual geometry, nonlinear material
sources and finite surface laws remain authoritative. This is the original
isothermal, equal-pressure vapor-equilibrium bridge with first-order endpoint
coupling. It does not infer immersed-liquid permeability, solute permeation,
pressure-driven osmosis, sorption heat or a temporal error guarantee.
Seven authored G1/G3/G4 regression groups specify distinct ideal/nonideal baths,
zero-net-uptake through-flow, identical-bath symmetry, independent finite phases,
sealed kinetic evolution, native geometry, depletion, work and cancellation
refusals. Changed Rust sources passed parsing; these integration regressions
remain unexecuted because the dependency build is blocked before fs-couple.

## Adaptive diffusion between two finite solution baths

`FiniteMoistureSolution::advance_paired_adaptive_to`,
`advance_paired_cartesian_adaptive_to` and `advance_paired_radial_adaptive_to`
complete a requested physical duration through the original first-order paired
endpoint and original full-versus-two-half controller. Every cell's moisture,
each bath's solvent mass and each bath's activity enter independent temporal
comparisons. Optional per-bath kinetic tolerances additionally compare every
original component's actual dissolved AND crystalline stocks; the original
phase identity/order and rate laws remain authoritative. Radial transport also
retains its mass-weighted shell-discrepancy norm.

Only accepted native half-step face transfers enter the complete-window water
accounts. Each bath closes against its own assigned faces even when opposite
exchanges leave almost no net specimen uptake. Final original-to-endpoint
component audits, separate bath balances, field/face conservation and an
independent compensated sum of represented bath/cell endpoint/start water
inventories must pass before any state publishes. A failed duration exposes no
accepted prefix. Last-step source receipts describe the final physical half
step; they are not a retained history of rejected probes.

`PairedSolutionMoistureTimeControls::maximum_field_evaluations` bounds actual
nested field proposals cumulatively across coarse, fine, rejected and failed
temporal trials. It supplements both endpoint-call and per-endpoint root caps.
The shared charge occurs before native surface/field admission; source and
phase owners retain their original work bounds. Supplied solvent brackets are
neither widened nor recentered, and source/domain refusals remain fatal.

This is a local asymptotic first-order error estimate, not SDIRK2 or a global
accuracy guarantee. Temperature and face assignments remain prescribed; the
original equal-pressure vapor-equilibrium bridge supplies no liquid immersion
permeability, solute permeation, pressure osmosis or thermal feedback. Six
G1/G3/G4 regression groups specify an analytic time-refinement ladder,
independent continuous through-flow and dissolution references, native
Cartesian/curved geometry, cumulative-work atomicity and cancellation. All
eight changed Rust files passed the pinned formatter/parser check. These
integration regressions remain unexecuted because dependency compilation is
blocked before fs-couple.

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


`ReactivePiezoelectricTransient::step_with_acoustic_ports_candidate` also retains
the original capacitor charges and winding fluxes in this same spatial solve.
Signed outgoing-wave work enters the enclosing body/reactive-storage/heat/source
balance once; prescribed moving supports supply their complete reaction work.
The returned report retains every acoustic port response. Existing APIs use an
empty port list and preserve their original electrical histories and behavior.
G1/G3/G4 regressions in `tests/piezoelectric/spatial_transient_load/reactive/acoustic.rs`
compare independent scalar RLC/guide dynamics, stored energy, conductor heat,
wave work, moving supports, empty-port equivalence and transactional refusals.

`TetPiezoelectricTransient::step_with_acoustic_thermal_loads_candidate` couples
those guide waves to the original finite resistor and ceramic heat inventories,
including optional moving supports and retained intrinsic material memories.
Signed outward acoustic work crosses the assembly boundary and enters the
joint energy balance exactly once; only conductor and intrinsic material
losses are deposited as local heat. The thermal report retains each guide's
wave response. Existing source-domain, caloric-inventory, candidate ownership
and immutable-solve reuse rules remain in force. Ceramic tensors and guide
fluid coefficients stay fixed during this path; it does not add thermal
feedback to reactive circuit storage or to the separate temperature-dependent
tensor path. The acoustic thermal regressions in
`tests/piezoelectric/spatial_transient_load/thermal/relaxation/acoustic.rs`
compare independent motion, field, memory and resistor-temperature equations,
separate finite heat cells, guide work, moving supports and atomic refusals.

### Finite-temperature full-tensor spatial piezoelectric dynamics

`spatial::transient::load::thermal::tensor::ThermalTensorPiezoelectricBody` couples the full native 3-D spatial fields and consistent-mass motion to one uniform finite ceramic temperature, with original per-cell tensor orientation and complete source-bound caloric inventory. Temperature proposals refresh the complete affine stress-charge C/e/epsilon source tensors and reassemble the same native operators under original budgets. The model may own refreshed tensors while borrowing unchanged geometry; identical-temperature trials reuse accepted operators. Each changed-temperature trial first solves electrostatics at fixed original displacement/free charge, restores physical momentum/history, and then runs the original shunt/current/support-motion midpoint solve. Signed reversible field/caloric transfer is derived from the full Helmholtz and entropy endpoints; actual resistor Joule heat is deposited only through its original receiver. Independent acceptance checks total field internal energy, kinetic energy, all finite caloric energy and actual external work/heat. Complete source/context/density/inventory checks, declared strain/electric-displacement bounds and a whole-interval caloric-plus-field capacity floor precede coupled Newton publication. This is fixed-reference geometry with one uniform ceramic temperature, fixed zero-voltage grounds and charge/current electrodes; intrinsic memory, nonzero ideal-voltage constraints, reactive circuit thermal feedback, phase transitions and independent cell temperatures require additional accounted models. Endpoint-temperature coefficients make the thermal feedback generally first order in time. Synthetic regression definitions cover moving/rotated supports, reversible caloric changes, heating-dependent capacitance at conserved charge across mesh refinements, separate finite shunt heating and transactional refusals; native execution remains a separate verification obligation.

`TetPiezoelectricModel::assemble_thermal_materials` retains the original resolved
thermal samples and their explicit expansion/clamped-pyroelectric sources.
Static reference forces and charges are assembled once from the actual weak
form. Static and transient fields recover `sigma=C*(S-alpha*dT)-e^T*E` and
`D=e*S+epsilon*E+p_clamped*dT` with the same reciprocal tangent. Reference
free energy at zero total strain/field is reported separately; ramp work
includes initial stress and removes initial polarization from transferred
charge. Fixed-temperature midpoint storage uses the shifted elastic potential
with retained-precision endpoint contractions. This reference datum is not
heat or an extra external energy supply. Zero-offset models retain the original
response path. Nonzero bias requires an explicit harmonic operating point.
Fixed-temperature transient memories retain the same thermal reference:
Maxwell strain increments are increments of S-alpha*dT, equal to delta_S at
held temperature; Debye polarization is additional to p_clamped*dT.

The finite thermal body admits the matching source/reference identities,
including at zero offsets, and uses these same bias-aware samples in every
changed-temperature trial. The caloric inventory represents stress-free,
zero-field storage; thermoelastic and electrocaloric transfers follow the
material potential's entropy and remain separate from dissipative heat.
Both total and shifted elastic strains obey the small-signal envelope.
`step_with_acoustic_ports_candidate` also receives/emits guide waves using
traces on the original model. Identical borrowed geometry/frame permits reuse
of the validated pressure/flow map across material snapshots without repeating
the whole-mesh face-owner scan. Signed wave work crosses the joint energy
boundary once, and the report retains the original acoustic responses.

Focused synthetic regressions cover open/short and free/clamped thermal
response, an independently solved shifted oscillator, reference/zero-bias
behavior, finite heat-induced voltage, signed electrocaloric exchange, and
secondary pyroelectric charge with moving-support/radiation work. Source,
amplitude and cancellation failures preserve accepted field and thermal state.

### Harmonic response about a steady piezoelectric bias

`TetPiezoelectricModel::harmonic_operating_point` solves and retains the actual
static boundary problem, including sourced thermal strain and polarization.
The immutable `HarmonicOperatingPoint` keeps exact model ownership, original
support/electrode topology and gauge-relative potentials. Its `solve_harmonic`
reuses the original sparse tangent/inertia/coupling and reciprocal acoustic
ports. Quadrature recovery is explicitly incremental: DC stress, charge and
reference storage never become phasors. The body mean free energy is static
free energy plus incremental mean storage, with no sinusoidal factor on DC.
Total strain, shifted elastic strain and electric field are admitted over the
entire offset elliptical cycle using a scaled two-dimensional dual bound;
quadrature-only samples and a loose triangle bound are not acceptance oracles.

R/RL steady current and external capacitor charge follow retained DC voltage
differences and original branch sources. Required ideal DC electrode supplies,
DC conductor heat, external storage and their AC averages are separate reports.
An ideal inductor requires an explicit persistent current and zero DC voltage;
its flux is not silently reset. This API holds the caller's DC bias rather than
solving an autonomous bias circuit. Thermal state and reference geometry remain
fixed; geometric stress stiffening, nonlinear DC tuning and temperature
modulation are outside the affine model.

`harmonic_operating_point_with_relaxation` solves the actual relaxed DC
operator with epsilon_infinity plus all directed Debye increments, while
Maxwell branch stresses vanish. It borrows the original immutable intrinsic
laws and reuses them for AC dispersion, memory and dissipation. Physical pole
charge is independently integrated into native field/electrode recovery;
static memory energy and polarization ramp work extend the same bias audit.
The DC memory subaccount is already inside total bias storage and receives no
sinusoidal factor. Zero-mean AC heat does not acquire a spurious DC loss.
Both zero frequency and each AC probe must lie in the original source bands;
a finite upper bound is admitted without authorizing unrestricted transients.
No new material card, tensor identity, transient history or second physics
trajectory is manufactured to represent this equilibrium.
Focused analytic regressions cover hot open/short response, bias-plus-cycle
limits, acoustic work, gauge-independent shunt work/storage and refusal.

### Settled-polarization spatial transients

`TetPiezoelectricTransient::new_with_equilibrium_relaxation` starts the
original midpoint body with fully relaxed Maxwell branches and settled Debye
polarization. One shared relaxed electrostatic solve retains the prescribed
initial displacement/velocity and original mechanical boundary topology.
Charge-controlled admission uses the actual relaxed field, with no intermediate
instantaneous-permittivity solution. Exact source/frame and broadband claims
and the complete cell scalar budget are checked before solving or allocating
histories. Native conductor charge and memory free energy are recovered once;
initial cell reports expose the retained inventory with zero heat/change.
The ordinary `with_relaxation` path still declares zero initial memory,
including under a fixed thermal preload.

Subsequent motion, moving-support work, current/shunt/reaction loading and
acoustic work retain the same coupled solve and independent energy accounts.
Total and shifted elastic strains remain inside the requested envelope.
A finite caloric receiver cannot change the temperature of a retained thermal
sample while its expansion and clamped-pyroelectric offsets stay frozen,
including a sample initially at its reference temperature. Such evolution
needs an accounted thermal memory owner; this increment is fixed-temperature.
Focused regressions exercise held hot bias, open-circuit voltage, relaxed-field
admission, polarized resistor discharge and Maxwell motion under thermal preload.

`release_electrodes_candidate` on the native and reactive transient owners
disconnects selected ideal voltage sources at an accepted endpoint. Each
released electrode retains its actual conductor charge; displacement, velocity,
time, fields, intrinsic memories and storage are unchanged. The reactive owner
also retains all attached source branches, capacitor charges and winding fluxes.
The native boundary map keeps current moving-support coordinates and original
electrode topology. At least one voltage reference remains per component.
Reference changes use retained relative potentials, with representability
admission before publication; large absolute gauges never enter that subtraction.
Subsequent steps use charge/current control at released electrodes. This
zero-time operation supplies no charge impulse or work and performs no new
field solve. It does not model a voltage-clamp jump, switch resistance or arcing.
Focused physical regressions cover partially charged Debye redistribution,
weak fields under large gauges, moving supports, reactive storage and refusal.

## Moving stress histories on the pressure/PCM path

`pcm_wav::observation::quadrupole::moving::MovingQuadrupoleRenderer` reads
original medium-frame Q/Qdot/Qddot at each source's retarded time and sums
native moving quadrupole pressure before audio conversion. Its existing
`PressureRenderer` interface feeds the original decimator and WAV writer.
Complete-window admission covers source support, continuous closest approach,
maximum Doppler clock rate, finite work and scratch budgets. Received
bandwidth must cover Doppler and geometric amplitude modulation below Nyquist;
a declaration does not establish source interpolation or spectral accuracy.

Callbacks allocate nothing, preserve global integer sample clocks, use DD
source summation and publish output/cursor only after complete success.
Out-of-support queries, failed reads and cancellation preserve accepted state.
Scope includes uniform and admitted quintic source translations, fixed or
moving microphones and the native complete ideal-plane/wedge image groups,
in the original homogeneous lossless fluid frame and f64 clock domain.
No source evolution, interpolation, moving quadrupole particle velocity,
directional microphone law or turbulence-model validation is inferred.
`tests/moving_quadrupole_pcm.rs` covers history/bandwidth/clearance refusal,
coherent cancellation, block invariance, atomic failure and native WAV parity.

`QuadrupoleMotion` retains each existing `CompactMotion::Observer` rather
than substituting a fixed microphone. The native direct and reflected owners
supply every path's original spatial Green Hessian, with source acceleration
and analytic jerk. Receiver motion affects the observation point and admitted
clock rate; it does not become source-clock stress derivatives. Readers keep
the physical source index and fixed medium-frame tensor, queried separately
at every image's own retarded time. Native inverse-image geometry and parity
apply the full rank-two reflection once. The raw DD tensor contraction is
summed over all physical sources and images before output rounding.

Admission counts complete source/image paths before retarded queries, charges
path-times-sample work, and includes each path in the lifetime root cap.
Every direct and reflected arrival must fit the declared source support.
The original whole-window Doppler bound covers the complete image group and
independently admitted source/receiver velocities. No callback is queried
during admission. `source_count` and `source_path_count` distinguish the
physical inventory and propagation work; `observer` retains original identity.
Per-event work is O(total paths) times the existing bounded root solve, with
one reusable block buffer and no extra trajectory/image bank.

Five G1/G3/G4 `tests/moving_reflected_quadrupole.rs` regressions cover moving
microphones with independently differentiated image-source potentials,
fixed/moving rigid and release wedges, preserved trace/shear and weak
destructive-interference residuals, full image support/work/Doppler admission,
and late image-read failure/cancellation with chunk-invariant retry.
All five, ten prior native moving/accelerating quadrupole tests and four
prior renderer tests passed focused native builds of the actual source
modules and original math/material/interval owners. The previous WAV-end-to-end
test and full Cargo/DSR workspace gates were not executed in that focused lane.

## Energy transport of original acoustic field observations

Moving monopole, force and quadrupole field renderers expose
`energy_transport(&field)` through the shared acoustic energy reduction. The
original admitted material and motion frame supply rho, c and laboratory mean
fluid velocity U. Source motion and microphone motion never substitute for U.
The input is one coherent physical p/u sample from that renderer, before
microphone projection or audio filtering; the query advances no history or
sample cursor and allocates no memory.

The quadratic perturbation energy is `E=(p^2/(rho*c^2)+rho*|u|^2)/2` and its
laboratory transport is `p*u+U*E`. Existing field samples retain their signed
pressure-work flux `p*u`, layout and callback behavior. The energy includes
all retained perturbation velocity, including a convected vortical wake; it
does not represent a change in background-flow energy or a radiated-power
qualification. This balance assumes the original homogeneous lossless fluid
and uniform background flow. Quadratic products can have twice the field
bandwidth, which matters when sampling their spectra or time averages.

## Moving stress fields and directional PCM

`quadrupole::moving::field::MovingQuadrupoleFieldRenderer` binds the original
`PreparedQuadrupoleHistory` to the admitted uniform/curved source motion,
moving receiver and complete ideal-image paths. Its explicit zero-prehistory
constructor requires authored initial stress, stress rate and stress integral
to vanish, and each trajectory to start at its original recorded position.
The original spline supplies Q/Qdot/Qddot without a second reconstruction.

The full velocity includes the spatially weighted past-worldline stress
integral required by linear momentum, not a pressure/impedance approximation.
It reuses the existing bounded vector quadrature and continuous past-clearance
admission. Fixed receivers cache complete source/geometry-cell prefixes;
moving receivers integrate at their current virtual receiver position.
Original physical tensors and history integrals reach the field owner before
it applies image parity and the polar-vector transformation exactly once.

Coherent pressure and Cartesian velocity are summed in DD before signed p*u
intensity and first-order microphone pickup. The same sample cursor supplies
full fields and the existing pressure/decimation/WAV interface. Callbacks
allocate nothing and publish output and cursor only after complete success.
History memory, original path/sample work and lifetime quadrature work have
explicit caps; failed and cancelled attempts do not refund numerical work.
Adaptive quadrature reports an error estimate, not an interval certificate.
The source compactness, bandwidth, spline convergence and prescribed-motion
validity remain producer obligations. Scope is the homogeneous lossless fluid
rest frame with fixed-axis ideal microphones and native ideal boundaries;
unknown initial fields, arbitrary walls and capsule dynamics are not inferred.

Five focused `tests/moving_quadrupole_fields.rs` regression definitions cover
stationary near fields and shear pressure nulls, block/cancellation behavior,
prehistory and finite budgets, moving receivers and oblique ideal reflections.
Native execution remains a separate verification obligation.

## Weak directional modal damping

`ModalPlateViscosity::from_psd_modal_matrix` retains positive dissipative
directions resolved by its compensated arithmetic even beside a much stronger
direction and an exact null. Its pivoted factor uses the existing double-double
owner for LDL Schur updates, with locally propagated arithmetic uncertainty.
Square roots and the arbitrary reference resistance enter only the realized
channels. Larger coupling ratios are formed before smaller cross terms to
avoid losing representable updates across widely separated scales. The former
`1e-10 * max_diagonal` positive-loss cutoff is absent; negative-pivot admission
and full reconstruction checks retain their existing boundaries. The ordinary
strictly positive-definite Cholesky path is unchanged.

This matrix-only API cannot identify positive source/projection error as an
exact physical null. Producers remain responsible for preserving their known
null spaces, and numerical Ritz vectors are checked by their actual physical
power rather than an assumed exact channel count. Focused native regressions
retain `1e-14` and `1e-100` loss ratios, correlated weak directions, dependent
rows, dense basis covariance, extreme-scale null motion and finite decay/heat
over a weak channel's physical relaxation time. These checks are ordinary
floating-point validation, not an interval or source-applicability certificate.

### Original stiffness in hereditary factorization

The shared regional relaxation factor supplies its original dimensional
stiffness entries to the compensated PSD fallback. Diagonal normalization is
confined to the established strictly positive-definite Cholesky path; a failed
Cholesky must not make rounded normalized entries the material source. This
preserves resolved weak strain energy and exact dependencies for regional
plate memory, shell substrate relaxation and coating relaxation, with the
existing full-form reconstruction and cancellation checks. Successful
positive-definite channel arithmetic remains unchanged.

## CST/DKT shell material damping

`ModalPlateViscosity::from_cst_dkt_shell` consumes the original CST/DKT shell
mesh, sparse mass pencil, reduced basis and one `ShellViscosity` per facet.
The spatial owner projects the passive membrane/bending resultant tensor,
including signed membrane-bending coupling, before the existing common owner
factors every retained modal cross term into dashpot channels. Physical strain
rates exclude numerical drilling stabilization. The supplied mass-normalized
basis is checked against the complete original mass matrix, including its
consistent and translation/director terms; no mode, frequency or normalization
is changed. Work admission covers projection, mass products, Gram checks and
factorization, with cancellation before publication.

The same channels enter `CoupledModalSystem` and its separate physical-loss
accounts. They introduce no spring, numerical damping, thermal destination or
new integrator. Source validity, frequency-band applicability, modal
completeness and use of the same original mesh remain caller obligations.
MITC and tied-shell kinematics retain their existing owners. The focused tests
compare sparse spatial power, the full modal quadratic form and channel power;
exercise lumped and consistent inertia, basis permutations, zero loss and
refusals; and compare coupled decay with an independently refined RK4 solution
of the original full modal equation while checking physical heat closure.

### Source-bound homogeneous shell viscosity

`plate_viscosity::shell::material::ShellMaterialViscosity` retains an immutable
`PlaneStressViscosity`, its original source state and plane, physical thickness
and orientation in the original facet. The existing material tensor owner
rotates the complete engineering-shear law and integrates homogeneous
resultants A = h V and D = h³ V / 12; B remains exactly zero. The preparation
refuses unrepresentable resultants and does not alter stiffness or inertia.
Material records are shared through `Arc` when multiple facets use one source.

`ModalShellMaterialViscosity::project` checks every retained eigenpair's
frequency against the intersection of all facet source bands, including zero
loss sources, then uses the original mass-normalized full CST/DKT projection.
The result retains facet sources, frequencies and full modal channels for the
existing physical-loss runtime. Cancellation and work refusal publish no
partial image. The source query point alone never establishes applicability;
bare mesh/model identity, homogeneous thickness correspondence, measured source
authority and modal completeness remain caller obligations. This adapter does
not infer layered, three-dimensional or hereditary material laws.

Focused tests compare rotated engineering-strain power independently, verify
linear/cubic thickness scaling, retain original records through modal
projection, and reject unsupported actual frequencies even for exact zero
viscosity. These are synthetic constitutive checks, not a material calibration.

## Material tensile admission on nonlinear acoustic strings

`NonlinearMaterialString::with_tensile_limit` explicitly binds one named static
strength claim per original homogeneous or perfectly bonded concentric layer.
The shared `AxialTensileCapacity` retains every source receipt, modulus and
annular geometry. No default strength, safety factor or allowable is inferred.
Missing, nonstatic or invalid strength refuses; measured whole-string and
helical winding sources require an explicit constituent stress model. Prony
sections and homogeneous hereditary bending sources refuse until their actual
layer memory stresses are included.

The existing Kirchhoff--Carrier storage supplies the actual common tension,
including geometric stretching from both transverse polarizations. Its exact
modal normalization supplies a continuous-span curvature envelope
`K = sum_n k_n^2 sqrt(2/(mu L)) hypot(q_xn, q_yn)`. Each layer is assessed using
`sigma_i <= E_i (T/EA + r_outer_i K)`, with the first maximum-utilization layer
governing. A single sine mode attains its curvature envelope; multimode
cancellation can make the envelope conservative. The two-plane norm keeps
this criterion independent of the transverse coordinate frame. The quantity
is a floating-point elastic model estimate, not an interval certificate.

Selected limits enter the original initial-state and candidate acceptance
checks. A crossed envelope refuses before physical state, material histories,
energy accounts or accepted time publish. Forces and stored energy remain the
original laws. Successful checks allocate nothing; preparation retains the
bounded modal factors and original capacity. Native reference, prepared and
pressure-rendering paths share the same check. Transversely moving supported
strings retain the guard through their original kinetic transformation: the
affine endpoint adds its actual stretching but has zero curvature. The
longitudinal moving-support extension below includes the original owner's
nonuniform axial stress. Transfer of a guarded potential
to the bow bank's geometric-only storage also refuses: that owner can add
independent axial material histories. Independently constructed bow banks do
not acquire this opt-in through source-property presence alone.

`tensile_assessment` exposes the selected elastic demand without enforcing its
limit, allowing diagnosis of overlimit candidates. This is a small-strain
normal-stress criterion at fixed material conditions and accepted temporal
endpoints. It supplies no viscous stress, torsional/shear failure, rate effects,
necking, fatigue, fracture energy, rupture topology or guarantee between time
samples. G0/G1/G3 regressions in `tests/material_string_nonlinear/tensile.rs`
specify axial-plus-bending demand, per-layer load sharing, transverse-frame
invariance, source refusals, unchanged admitted mechanics and native rollback,
including supported conversion. Native execution is a separate verification
obligation.

### Tensile admission on material fiber ribbons

`TensionedFiberRibbon::with_tensile_limit` binds every original fiber layer's
named static strength to the existing nonlinear `FiberStrain` owner. Its
retained longitudinal strain supplies the axial upper bound, including the
maintained preload and geometric mean. Original mass-normalized transverse
coordinates supply a continuous-span curvature envelope. Elastic layers use
`E_i (epsilon_max + r_outer_i K)` with their own physical fiber radii: refining
the ribbon width changes neither stress nor utilization.

Sourced Kelvin--Voigt layers retain independent original axial and bending
viscosities through `KelvinVoigtTensileCapacity`. The original longitudinal
owner supplies the complete strain rate, including `y_x*y_xt` and its mean.
For each layer, signed `E_i*epsilon_j + eta_ai*epsilon_dot_j` axial harmonics
and `E_i*kappa_n + eta_bi*kappa_dot_n` bending modes combine before their
span magnitudes are summed. Unloading can therefore cancel tensile stress;
no section-average viscosity redistributes it between layers. The maintained
preload has zero rate. `tensile_state_assessment` reports total stress, while
`tensile_assessment` explicitly reports the elastic contribution alone.
Both retain the original section/strength receipts and force per physical fiber.

Preload and the original coupled candidate gate refuse a crossed envelope;
state, histories, work accounts and time remain transactional. Forces and energy
are unchanged. Velocity-dependent limits use the actual original modal velocities
at cold binding, component-state restoration, accepted states, converged dynamic
candidates, and final staged static equilibria, before publication. Cold restores
also admit all displacement-only constitutive domains and retain their accepted
material histories. The held-force displacement increment
divided by dt is not substituted for an endpoint velocity. Mapped admission
buffers are prepared once and reused; displacement-only laws gather no extra
velocity data. The fiber assessment uses bounded stack workspaces and cached
per-layer coefficients without new heap allocation or runtime source lookup.

This remains a fixed-state normal-stress envelope at temporal endpoints, with
no guarantee between samples, rate-dependent strength, yielding, fatigue or
rupture topology. Hereditary combinations refuse in either builder order until
actual layer memories can participate. Tests in
`tests/material_multi_bow/fiber_ribbon/tensile.rs` and its `tensile/viscous.rs`
child cover analytical demand, width refinement, distinct layer viscosities,
signed cancellation, source/preload refusal, unchanged admitted mechanics,
and native cold/restore/static/dynamic admission and rollback.


## Actual worn material leaving hardness-coupled thermal receivers

`PlasticSlidingHeatContact::export_wear_candidate` transfers newly accumulated
original A/B Archard volumes to the existing
`DepositionSimulation::remove_material_before_pending_events` owner. Each
receiver supplies its original sourced caloric density at unchanged current
specific enthalpy. Native inventory checks require that density to agree with
stored mass and volume; temperature-dependent density is not silently replaced
with m/V and no thermal-expansion law is inferred. The exact caloric objects,
original material/body IDs, pressure/indentation/shear/wear histories, and clock
survive. Removed mass carries its current signed enthalpy, including negative
datums; remaining specific enthalpy and temperature do not change.

The caller authors complete remaining prism geometry, every incident contact
patch, and current original conductivities. Admission retains the same contact
link, ordered cells, actual loaded area, original prismatic face areas and axis
registration. Only normal shortening/translation is supported; rotation,
tangential translation and a silently enlarged receiver face refuse. The next
coupled heat interval uses actual reduced masses and new half-cell distances.
Original cumulative wear-fraction and homogeneous-depth envelopes are retained.

A per-body checkpoint records the cumulative requested wear already represented
by the supplied thermal binding or an accepted export. Pre-binding wear is never
removed retroactively; repeated export requests zero additional material, and
later sliding exports only its new wear. The native report separately retains
requested versus represented volume and mass discrepancies, outward enthalpy,
and stable inventory/energy balances. The contact and thermal successors must
be adopted together. Any source, geometry, density, patch, representability or
cancellation refusal leaves the inputs untouched; pending endpoint events stay
pending.

This is an explicit heat-then-export composition requiring time refinement,
with caller-supplied uniform thermal control volumes. It does not reconstruct
a spatial thermal mesh from local wear, model debris evolution or surface
creation energy, or account for mechanical work needed to move those volumes.
G0/G1/G3 regression definitions cover separate material inventories and
enthalpy, an independent two-capacity heat oracle after export, unchanged
mechanical history, no duplicate or retroactive removal, original cumulative
bounds, and transactional patch/registration/cancellation refusals.

### Adaptive prescribed sliding with hardness, wear and finite heat

`PlasticSlidingHeatContact::march_prescribed_velocity_candidate` refines a
constant actual per-cell velocity span with the original normal boundary held.
A full interval and two native half intervals start from the same accepted
history; only the half-step candidates and their original source payloads may
advance. Componentwise differences control thermal-cell energy, interface
temperature, pressure, approach, indentation, each body's recession, accumulated
slip and physical work. Relative scales use changes from the interval origin.
No extrapolated material, contact or thermal state is manufactured.

The adapter uses zero tangential body compliance. Subdividing a released body
displacement while retaining nonzero compliance would change the motion problem
and is not supported. Whole-span contact, thermal and combined work gates are
independent of per-step residual tolerances and substep count. Compensated sums
use native incremental storage and preserve signed parameter work, friction,
direct plastic heat and nonthermal plastic work separately. They do not supply
a full body thermodynamic or passivity theorem.

Attempt and accepted-receipt budgets bound work and storage. Source, numerical,
coverage, inventory or cancellation refusal publishes no partial span, including
after private accepted intervals. Native failures remain explicit; only excessive
coarse/fine defects trigger refinement. No thermal event or material export is
consumed, and cumulative removal limits/checkpoints survive. The defect is a
local estimate for the first-order coupled split, not a global error bound.

### Spatial magnetic hysteresis and finite enthalpy

The opt-in `hysteretic-eddy-heating` feature connects the original
`PlanarHystereticEddy` vector-history owner to `MagneticThermalNetwork` and
the existing finite-cell enthalpy Newton driver. Every thermal trial starts
from the same accepted potential and vector polarizations. Native per-triangle
pinning, intrinsic viscous switching and spatial Joule heat enter their actual
receiving cells exactly once. Insulating magnetic triangles can heat through
pinning and finite-rate switching without Ohmic conduction; any
dissipating triangle without a finite receiver refuses. The endpoint field's
numerical dissipation remains in the work account and never heats a cell.

Actual endpoint temperature refreshes the original conductivity source.
Every mapped magnetic fit must be source-constant over the complete original
caloric temperature interval, with the same material card, selected claims,
quantity schema and fixed external conditions. The report retains complete
fitted magnetic receipts separately from the electrical bundles. Constant
solid density and whole-cell volume correspondence remain mandatory. Stable
magnetic storage change, physical heat, numerical loss, thermal storage and
source work share one acceptance audit; failure/cancellation publishes neither
history. The shared report's `deposited_joule_j` counter denotes total deposited
physical heat in this driver; `deposited_heat_j()` names that sum explicitly,
and its native field report separates Joule, pinning and viscous heat.

An explicitly source-selected viscosity uses the same physical interval as the
thermal driver and every native field trial. The finite receiver deposits the
full constitutive heat, never the half-viscous variational objective. Its source
bundle includes all original viscosity claims, even an all-zero profile, so the
complete caloric interval must support the kinetic fit too. Native winding
circuits share this existing eddy and thermal route. Synthetic regressions
check finite-rate polarization and continued hold heating against an independent
scalar kinetic balance, exact deposited energy, retained source queries and
atomic refusal of insufficient viscosity temperature support.

This is fixed-geometry, fixed-fit, first-order magnetoquasistatic/thermal
coupling. It supplies no changing-fit thermomagnetic entropy/parameter work,
moving media, phase change, fitted core-loss overlay or experimental material
validation. Focused synthetic regressions compare insulating minor-loop
heating and conducting temperature feedback with independent scalar play and
P1 balances, distinguish numerical loss from heat, and exercise missing
receivers, source-temperature support, volume, cancellation and work budgets.

### Retained dielectric switching in native circuits

The separate opt-in `experimental-dielectric-hysteresis` feature exposes
`material_circuit::dielectric_hysteretic::DielectricHystereticCircuit`.
One original capacitor is replaced by a private native voltage port. Every
charge trial retains the same accepted source-bound polarization memory;
the exact monotone-segment mean voltage closes `delta q = dt I` through the
original Kirchhoff DAE. Native R/L/C coordinates, source ports and material
resistor bindings retain their owners. Private work cancels internally, while
resistor heat and dielectric switching heat remain distinct and enter the
combined energy audit once. Failed trials never publish charge or polarization.

Three instantaneous-current probes refuse detected charge reversals. They do
not certify that an arbitrarily large interval contains no hidden turn; the
caller supplies waveform/reversal cuts and refines threshold crossings. The
zero-branch limit is the native midpoint circuit; fixed active-set RC charging
has an independent analytic refinement regression. Geometry and coefficients
remain fixed. Reported switching heat is not automatically deposited into a
thermal cell, and this owner supplies no leakage, electromechanical strain,
spatial electric solution or empirical calibration. Six native integration
regressions cover linear RC/RLC limits, exact loop work, charging refinement,
atomic refusal, detected LC turning points and small transfers beneath large
unrelated energy biases. Both native and combined audits use the original
solver's stable energy increment.

### Finite thermal feedback for retained dielectric switching

With the same default-off `experimental-dielectric-hysteresis` feature,
`material_circuit::dielectric_hysteretic::thermal::ThermalDielectricHystereticCircuit`
couples accepted charge/polarization to the existing finite-cell enthalpy driver.
The original `PreparedMaterialNetwork` resolves caloric storage, resistance and
conduction at each proposed temperature. `HeatAssembly` providers retain that
bulk network and can supply thermal contacts, layers, closed radiation and
finite coolant. Every trial restarts the SAME accepted electrical history;
this adds no thermal, constitutive or nonlinear solver.

Switching work enters the actual whole fixed solid dielectric volume; native
resistor work enters its original material cell. An explicitly authored leakage
resistor may share that receiver: both physical losses accumulate once. Their
separate receipts survive, while `deposited_heat_j()` reports the total.
Numerical dissipation is never deposited as heat. Original circuit, finite-cell
thermal, total-energy and closed-radiation gates remain active.

Source-bound switching requires the exact original caloric card and accepted
query, plus explicit support for the entire admitted temperature interval in
every selected constant coefficient. Each trial refreshes that same fit; source
query, polarization, charge and thermal clock publish together. Both actual
cell volume and mass/density must reproduce capacitor area times gap, with a
wholly solid constant-density caloric chart. Cancellation, source domains and
combined nested circuit-work budgets refuse the complete proposal.

Endpoint resistance with backward-Euler heat is generally first order; charge
reversals and source/thermal events require interval cuts and refinement.
Electrical contact work and field-dependent thermal boundaries remain outside
this lane. No changing switching coefficients, electrocaloric heat, geometry
evolution, breakdown or empirical accuracy is inferred.
`tests/dielectric_hysteretic_thermal.rs` supplies independent two-cell
leakage/switching/R(T) balances through the native heat assembly, major-loop
finite heating, source/history continuation and atomic refusal regressions.

### Spatial dielectric hysteresis on the finite enthalpy timeline

The default-off `experimental-dielectric-hysteresis` feature also enables the
native field lane and exposes
`circuit_enthalpy::dielectric::hysteretic::ThermalHystereticElectrostatics`.
It composes the original planar Gauss/electrode history owner, source-bound
material network and finite-cell enthalpy driver; no constitutive, field or
thermal solver is duplicated. `DielectricFieldSolve` accepts the native single
endpoint controls or the existing adaptive affine-path controls through the
same thermal entrypoints. The complete original HeatAssembly may supply
conduction, thermal contacts, layers, closed radiation and finite coolant.

Every original triangle has an explicit `DielectricRegion::Cell` or `Fixed`
correspondence. Mapped triangle volumes, actual receiver volume and mass/density
must agree within the declared absolute geometry allowance; multiple triangles
may partition one whole receiver. Its source must be the same original caloric
card and complete accepted query. Fixed geometry requires a wholly solid,
constant-density caloric chart, and every selected constant switching
coefficient must explicitly support that chart's entire temperature interval.
Any nonzero physical switching heat in an unmapped Fixed region refuses.

Fixed fits and prescribed loading make the electrical endpoint/path independent
of temperature. It is evaluated once before the thermal residual, with the
native field/path work budgets retained. Each thermal trial rechecks source
queries and transport while borrowing the SAME immutable electrical account;
it neither reruns Gauss nor copies or advances polarization memory. Per-cell
physical switching heat is accumulated once into receiving-cell powers.
Numerical electrical dissipation is retained separately and is never heat.
The original field/source-work, finite thermal, total-energy and closed-radiation
gates precede publication of both histories, final source queries and clock.
Cancellation, work budgets and source-domain refusals preserve accepted inputs.

The optional path summary retains every coarse/half/rejected field solve and
operator application, plus accepted-only empirical refinement estimates. Both
the candidate getter and consuming parts preserve it. Accumulated electrical
path heat Q is deposited as Q/dt over ONE backward-Euler thermal interval;
electrical refinement does not resolve heat-release chronology or increase
thermal time order. Supply waveform/event cuts and refine the thermal interval.
There is no changing switching fit, electrocaloric transfer, field-dependent
thermal boundary law, leakage or geometry evolution in this lane.

`tests/hysteretic_dielectric_enthalpy.rs` supplies six focused G1/G3/G4 fixtures:
independent scalar-play heating in a rotated multi-triangle receiver; independent
two-cell temperature-dependent conduction through HeatAssembly; Fixed-region
zero-loss/refusal behavior; original source and whole-volume admission; late
cancellation, budgets, caloric support and replay; and cumulative refined-path
heat with separate numerical loss and retained global-work receipts.

### Resolving dielectric heat-release chronology

The same default-off spatial dielectric owner exposes
`advance_affine_candidate` and `advance_affine_candidate_with_heat_assembly`
through its `trajectory` module. One authored affine voltage/charge and volume
charge segment is now refined in physical time. Full and two-half trials use
the original coupled endpoint solver; only the two-half field and thermal
histories can be retained. Each substep deposits its own physical switching
heat while the original conduction, reservoirs and complete HeatAssembly act.
A hold can therefore refine solely for thermal cooling with unchanged electric
field, zero switching heat and unchanged polarization.

`HystereticThermalPathControls` combines the original P/Q/work refinement
policy with a separate L1 cell-energy defect and explicit global thermal-driver
evaluation cap. Thermal relative errors scale interval energy changes, never
absolute reference enthalpy. Both per-cell heat defects and thermal defects
use absolute sums so opposing spatial errors cannot cancel. Coarse, half and
rejected trials consume one global field-solve/operator/evaluation budget.
The original per-step caps remain active. Duration-weighted native energy
allowances and thermal Newton absolute targets precede whole-path physical
balance gates; no converged-state or small-update fallback bypasses them.

Only interval descriptors and a bounded set of trial states are retained. A
rejected first half is reused only for its identical left-child starting
history. Cumulative accepted-only electrical/cell and thermal energy accounts
use compensated sums; the original stable thermal residuals are preserved.
The report explicitly identifies its source sample, nonlinear reports and
charge residuals as final-substep data. Global attempted-work counts and
empirical error estimates remain in `HystereticThermalPathSummary`.
Cancellation, refinement/work refusal or source-domain failure returns neither
new history. These are first-order thermal trajectories with empirical local
estimates, not certified global time errors, automatic event detection,
temperature-dependent switching or an increase in the thermal integrator order.

The three `hysteretic_dielectric_enthalpy::trajectory` G1/G3/G4 tests compare
switching during cooling with an independent piecewise analytical solution,
check exponential cooling through HeatAssembly while every electrical error is
zero, and exercise global work caps, depth, late cancellation and exact retry.


## Accelerating stress histories through the original PCM renderer

`MovingQuadrupoleSourceWindow<M>` and `MovingQuadrupoleRenderer<F,M>` now
accept uniform motion, a borrowed `PreparedSourceTrajectory`, or the existing
mixed `CompactSourceMotion` ensemble. `QuadrupoleMotion` extends the sealed
original motion owners with complete stress-potential derivatives; the
renderer retains one source-history, admission, coherent summation and atomic
streaming loop. Its default remains `UniformSourceMotion`.

Curved paths use their original quintic kinematics and analytic jerk at each
retarded source event. Their continuous clearance and window-local Doppler
bounds come from the existing `CompactMotion` implementation. Observer, medium,
limits, source support and receiver coverage must match. Each original root
keeps its finite refinement cap; `maximum_retarded_iterations()` reports the
checked complete-window cap. Admission and constant-work endpoint/final checks
are additional bounded work. Runtime trajectory lookup and the prepared
18-coefficient pressure contraction allocate nothing in the renderer.

The receiver remains fixed and omnidirectional in a homogeneous quiescent
lossless medium. The declared receive bandwidth must account separately for
geometric modulation: the Doppler clock-rate bound does not certify a curved
source's full emitted spectrum, especially across C2 knots with jerk jumps.
History reconstruction, trajectory convergence and antialias accuracy remain
producer responsibilities. The existing decimator/WAV path consumes the actual
coherent pressure without a second pitch transform or receiver multiplier.

G3/G4 `tests/accelerating_quadrupole_pcm.rs` specifies mixed-owner retarded
clocks and chunk invariance, observer/history/Doppler/work refusal, late failure
and cancellation atomicity, and parity through the native decimator and WAV
writer. Native execution remains a separate validation obligation.

## Reed channel slit law (`slit_length`)

The optional reed-card row `slit_length LENGTH_M` selects `ReedBoreVoice::set_viscous_slit_length`, which builds `ViscousSlit::separated_jet_from_gas`: the separated-jet entrance law of van Zon, Hirschberg, Gilbert and Wijnands 1990 (J. Phys. Colloques 51 C2-821, HAL jpa-00230506, eqs. 1-8 and Fig. 2) on the voice's exact inlet gas. `LENGTH_M` is the reed channel's streamwise length (Chatziioannou and van Walstijn 2008 print 1.4 mm for their mouthpiece), not the lay length of the reed. With `v = phi / (w h)`, `Re = v h / nu` and `S = 2 h` (the paper's fit, within 3 % of its 2-D measurements for `Re <= 4000`): below `L / h = 3` the jet never reattaches and `dp = rho v|v| / (2 alpha^2)`, `alpha = 0.45` (eq. 3); otherwise `dp = rho v^2 / (2 (1 - d)^2)` with `f(d) = 4 d + 9 ln(1 - d) + 5 d / (1 - d) = 6 (L - S) / (h Re)` while the boundary layers are still growing (eqs. 4, 5, `d <= d_c = (4/9)(1 - sqrt(5/32))`), and the Poiseuille tail `dp = 0.743847 rho v|v| + 12 mu (L - S) v / h^2` past `l_c` (`l_c - S = c Re h`, `c = f(d_c) / 6`). Three print checks are in the code doc: eq. 4's braces are positive (the `5 b / (b - h)` of the 2008 paper is the misprint), eq. 8 lacks its 9 (`c` is taken from eq. 4 at `d_c`), and eq. 7 as printed is inconsistent with eq. 5 at `l_c` (0.903 against 0.731), so the tail is a construction from eqs. 1, 5 and 6 that reproduces the Fig. 2 solid curve at fifteen read points to 0.011 (test tolerance 0.015). The law is quasisteady, laminar and 2-D; the step at `L / h = 3` is the paper's two data sets with no blend; it refuses `linear_gap_flow` and `fixed_opening_pressure_drop` (uniform reed gap only). Every other constructor (`ViscousSlit::new`, `from_gas`) is still the fully developed law `dp = 12 mu L v / h^2 + rho v|v| / 2`, so the spatial, modal and material-boundary users do not move, and a card without `slit_length` builds no slit at all.

Tests: `viscous/separated_jet.rs` (constants against the printed closed forms, the series against the closed form of the braces, the thickness solve over twelve decades, the Poiseuille limit, `alpha` below `L / h = 3`, value continuity at `d_c`, the 15-point Fig. 2 table, tangents against differences, the speed solver with a wave load and a staple recovery on all three branches) and `viscous/jet_tests.rs` (the public entry points against a bisection reference of the printed law, the characteristic junction against manufactured flows, the staple junction's equation and its `C_P -> 0` limit, the refusals) and the voice-level same-sample closure tests in `render/scattering.rs`.

## Reed characteristic-line ladder row (`line_ladder`)

The optional reed-card row `line_ladder TAIL_TOL MAX_TAPS`, after the hole rows (`open_hole_interaction`, `tonehole_lengths`, `constriction_noise`, `open_hole_orifice`) and before `compile_limits`, declares the convergence ladder of the reed's characteristic line (`RegionalDuct::with_line_ladder`; `LineLadder` with `peak_tol_cents` at its default 0.05). `TAIL_TOL` is the relative tail energy in the last quarter of the window (positive and below 1; default 1e-10) and `MAX_TAPS` the ceiling, a power of two in [256, 131072] (default 65536; the row's own bound, under the code's `MAX_LINE_TAPS` of 262144). A card without the row takes the default ladder bit for bit (`characteristic_line_response_with_poles` calls `characteristic_line_response_with_poles_ladder` with `None`; the tabulated-load route of `regional/baked.rs` follows the declared ladder too); a molecular bore refuses the row. A line that does not converge refuses as `AcousticRealizeError::ReedLadder { max_taps, tail_energy, tail_tol, peak_shift_cents, peak_tol_cents }`, whose text prints the measured figures of the last rung (it used to say only that the ceiling was reached). The row raises the ceiling or declares the tolerance; it cannot repair a response with a jump, whose tail stays near its value when the window doubles (the monopole TMMI load of the saxophone lattice card C4: tail energy 3.52e-4 at 65536 taps and 3.92e-4 at 131072, against 1e-10). A tolerance is a numerical choice: declare it before the render and state it on the page.

Tests: `the_line_ladder_row_is_validated_and_the_declared_default_changes_no_bit` (bad rows refuse; `line_ladder 1e-10 65536` renders the bits of a card without the row), `the_line_ladder_row_reaches_the_ladder_and_the_refusal_carries_its_figures` (`line_ladder 1e-300 512` refuses with `max_taps=512`, the measured tail energy, `(tol=1.000e-300)` and the peak shift).

**`open_hole_interaction tmmi 2 piston`**: the word `piston` after EPSILON, which must be 2 (any other value refuses), declares `fs_duct::OpenHoleInteraction::with_piston` (see the fs-duct contract: the monopole mutual impedance times the baffled-piston directivity of both mouths, the large-distance limit of the exact uniform-piston pair, not the exact integral). Without the word the card is unchanged to the bit, provenance included. The receipt names it only when declared, as `outlet_receiver.open_hole_interaction_mutual` (`ReedPerformance::open_hole_interaction_is_piston`), next to `open_hole_interaction_fallback`. Tests: `the_piston_word_is_validated_and_off_by_default`, `reed_outlet_record_names_the_piston_mutual_only_when_declared`.

## Double-reed staple pressure recovery and jet contraction (`staple_recovery`, `contraction`)

`ModalAperture` also composes a fully developed `ViscousSlit` with the same
explicit `StapleRecovery`, selected in either order before stepping. The full
mechanical response condenses to `D p_blade - R + Z q = D drive`, where
`D = 1 + Z pressure_mobility`. The channel keeps all of its linear no-slip
wall loss; recovery acts only on the quadratic jet head. Both the midpoint and
exact structural propagators, including a full material viscosity matrix, use
the same scalar solve and original pressure/volume work account. The blade
sees `p_blade`, the bore sees `p_blade - R`, and the signed exchange `-R` times
swept flow remains explicit. Forward jet loss is nonnegative in the original
recovery domain; this is not a new whole-system passivity claim. Zero recovery,
reverse flow and a closed gap retain the original viscous pressure arithmetic.
Spatial closures, independent vena contraction and the separated-entrance
channel law are not combined with this modal viscous recovery path. No new
material defaults or fitted coefficients are supplied. The focused tests in
`modal_aperture/{viscosity,propagator}.rs` compare manufactured full-matrix
momentum and work, the fixed-gap scalar limit, closure/reopening and the
zero-recovery trajectory; the scalar-law tests check condensation and refusals.

The DG-valve file's typed lumped `reed` now accepts the existing
`slit_length LENGTH_M` record. A positive streamwise length in metres selects
this fully developed laminar channel law, with viscosity and density from the
complete actual DG inlet gas, including an authored temperature/composition
field. It can compose with the reed's explicit `staple_recovery`; independent
`contraction` remains refused. `slit_length none` retains the original inviscid
reed. The existing finite-gas transaction refreshes the same source-bound
channel together with the acoustic gas; it preserves the authored length and
recovery and introduces no new thermal storage. This DG selection does not
select the separate reed-performance route's separated-jet entrance model.
The focused DG reed/card tests exercise source-gas binding, changed flow,
callback invariance, geometry/contraction refusals and finite-gas refresh.

Two optional reed-card records change how a massive reed's opening passes its flow. Both are off unless declared; absent, or at their identity value (`C_P` 0, `C_c` 1), the voice is the original one, bit for bit.

**`staple_recovery C_P S_C_M2`** (Almeida, Vergez and Caussé 2007, J. Acoust. Soc. Am. 121(1) 536-546, arXiv physics/0607011v1 pp. 8-9, eqs. 14-17), right after `slit_length`, or right after the `reed` row on a card without one, selects `ReedBoreVoice::set_staple_recovery`. A double reed's channel widens from the constriction of area `S_c` into the conical staple, and a fraction `C_P` of the dynamic head at the constriction comes back as pressure. The flow through the opening follows the slit law on the constriction drop `(dp)_c = p_m - p_c`: the viscous slit law with `slit_length`, the inviscid Bernoulli law without it (Almeida's own model neglects the viscous loss along the reed). The blade is driven by `(dp)_c` (eq. 16), not by the total mouth-to-bore drop (`JunctionPressure::blade_drop_pa` in `reed_bore::contact`), and the bore input is `p_r = p_m - (dp)_c + C_P (1/2) rho (q/S_c)^2` (eq. 17), `q` being the flow through the opening. `ViscousSlit::characteristic_junction_with_recovery` and `moving::characteristic_junction_with_recovery` solve the same-sample characteristic junction on the two laws (the inviscid one is the viscous one with a negligible slit, to 1e-9): for forward flow the quadratic term takes the factor `1 - C_P (S/S_c)^2`. Recovery acts on forward flow only (in reverse the staple is upstream of the jet) and on the JET flow: the moving blade's swept flow and any body flow are not carried through the staple. The law is singular where `C_P (S/S_c)^2 >= 1`, where it would cancel the whole jet loss; there it refuses (`STAPLE_RECOVERY_DOMAIN`), and `set_staple_recovery` refuses a reed whose rest opening already is there. The only printed numbers are one measured reed's, `C_P` 0.8 (Fig. 15) and `S_c` 4e-6 m^2 (Fig. 14): a ready-to-play natural-cane oboe reed, humid air, masses added, no lips. A card that uses them borrows them for another reed, and says so. Inside the domain the element dissipates `(dp)_total q`, never negative; because the blade feels `p_c` while its swept flow enters the bore at `p_r`, the step's energy account carries one exchange term, `-(C_P (1/2) rho (q/S_c)^2)` times the swept flow, of either sign.

**`contraction C_C`**, right after `staple_recovery` (or after `slit_length`, or the `reed` row), on a card without `slit_length` and with a massive reed, selects `ReedBoreVoice::set_jet_contraction`: the jet leaves the opening with the area `C_c w h`, so the inviscid flow is `q = C_c w h sqrt(2 (dp)_c / rho)` (the `rho U|U| / (2 (C_c w h)^2)` of `fs_phs::GapFlowLaw` without its channel terms) and the staple law takes `C_P (S/S_c)^2` on that jet area; the moving blade's swept flow is geometric and keeps the full width. `moving::characteristic_junction_with_recovery` takes the coefficient (1 is the original junction, bit for bit; outside `(0, 1]` it refuses). There is no default, and **no printed value for a double reed is held**: the clarinet's 0.6 (Hirschberg, in Silva, Kergomard and Vergez 2008) and the two-dimensional free-streamline bound pi/(pi+2) = 0.611 are borrowed numbers; Almeida's flow fit used no contraction, so an area taken from it may already hold whatever contraction there is, and a value given here is an upper bound on its effect. Hirschberg's downstream-neck resistance, which Vergez et al. (2003) carry to double reeds, moves the characteristic the other way from the measured oboe (Carral, Vergez and Nederveen 2011). The channel's unsteady inertance (`GapFlowLaw::inertance`) is not part of this record: no printed effective length of a double reed's channel is held.

Tests: the flow maximum of the printed characteristic sits where the closed form of eqs. 15-17 puts it; the domain refusal by name; the junction's power balance and passivity; the inviscid junction against the viscous one with a negligible slit; the contracted junction against the manufactured flows and against `GapFlowLaw::volume_flow` with the same contraction; the contraction moving the staple domain to the jet area; the same-sample flow closed on the contracted area; the blade momentum under `(dp)_c` and the step's energy account; bit-identical absence for both rows; callback invariance; the card rows' refusals.


## Source-bound plastic shell and spatial enthalpy feedback

`shell_plastic_enthalpy::ThermoplasticShell` binds each physical layer of every original CST/DKT facet to a unique active enthalpy cell, its original caloric object, reference volume and sourced density. The native plane-stress section retains independent complete material histories at all thickness points of all three surface points. Source families resolve actual thermal expansion, affine energetic modulus/hardening changes, yield softening and selected Voce/Perzyna/Chaboche kinetics. The native material account includes coefficient free-energy work and mechanical entropy; physical/reversible heat is integrated conservatively to the original cells and supplied exactly once to the existing conduction/observer timeline.

`advance` and `advance_backward_euler` share frozen-history equilibrium and thermal proposals. Every trial starts from the same accepted coordinates, complete plastic histories, enthalpies, observer state and, for dynamics, velocity. The actual interval owns both rate-dependent material response and native momentum. Only matching temperature/heat and independent force/support, thermal and combined first-law gates publish the shell, heat, motion and observer candidate together. Signed transverse residual work stays explicit. Material numerical loss, retained reference drilling storage/loss and physical/gauge kinetic integration loss never heat a cell. The original consistent mass is prepared only when requested, retained across intervals and checked against each original layer density, not merely total mass.

Physical sources and facet bindings are shared immutably; borrowed endpoint sections avoid material-law clones, and reference mass/gauge operators outlive temporary thermal laws. Source-bound conductivity selections lower to the existing material-network observer, retaining its events, reservoirs, conductances and interval hooks. The reduction remains fully solid, fixed geometry/density, small strain and rotation, with one temperature per physical layer/facet and explicit thickness quadrature. Finite forming, melting, transverse-shear plasticity, evolving transport geometry and delamination are not claimed. Physical interval, thermal substep, spatial cell and thickness refinement remain independent accuracy requirements.

Four focused synthetic regressions cover an independent heated uniaxial plastic membrane/permanent-set oracle, layer-gradient thermal moment reversal with original conduction, late-refusal checkpoint preservation, and prescribed plastic motion with an independent P1 kinetic integral, incorrect eccentric mass rejection and no heat from numerical stopping loss. They also exercise a cold reference binding whose free-coordinate budget is smaller than its total physical coordinate count. Formatting and source checks do not substitute for execution of these native numerical tests.

## Spatial thermoplastic indentation hardness

`solid_enthalpy::plastic::hardness` binds one independently selected original
indentation test to each tetrahedron of an accepted `ThermoplasticSolid`.
`bind_hardness_field` bounds point count and retains original mesh, ordered cell,
caloric source and test identity. `at_accepted` observes each actual caloric
temperature and native accumulated equivalent plastic strain, retaining complete
source receipts in original element order. It preserves full J2/Voce/viscous and
N-branch Chaboche histories, including reverse loading and thermal branch transport.
There is no averaged temperature, yield-to-hardness conversion, repeated return
map, repeated mechanical coefficient refresh or second mechanical history.

`check_candidate` admits an existing native thermoplastic proposal and returns
that unchanged proposal with its complete hardness field. Source coverage,
inventory, clock, original-mesh and cancellation refusals publish no partial
field or advanced physical owners; callers adopt both outputs together. These
are endpoint observations, with no assertion of coverage between physical
times. Consumers still own spatial projection, contact/indentation scale and
avoiding duplicate resolved and reduced indentation energy; no additional
pressure cap, wear law, annealing kinetics or heat term is introduced.

The focused synthetic `solid_plastic_enthalpy::hardness` regressions specify
heterogeneous element responses against an independent uniaxial plastic/heat
oracle, complete cyclic branch and native work preservation, and late endpoint
coverage refusal with accepted-owner preservation. They also exercise bounded
point count and original-mesh/clock admission. These are numerical source/history
checks, not measured material calibration.

## Ductile initiation on coupled thermal plastic solids and motion

`solid_enthalpy::plastic::ductile` adds `advance_with_ductility` to the original
`ThermoplasticSolid` and `ThermoplasticDynamics` owners. The existing native
proposal retains all thermal, equilibrium, source refresh, physical/reversible
heat, inertia and energy-balance behavior. Only after that complete proposal
converges does `DuctileHistory` assess the actual endpoint stress and accumulated
plastic alpha. Complete J2 and Chaboche histories remain the original states;
backstress is not projected away and reverse flow accumulates ductile usage.

Each point's physical temperature is evaluated from its retained cell enthalpy
through the SAME original caloric source. The shared physical duration is the
target time minus the accepted coupled time, including after restart. Material
kinetics, momentum and ductile plastic-strain rate use that interval. No ambient
temperature, iteration time, total strain or kinetic/numerical loss substitutes
for those quantities. The bridge allocates one bounded temperature vector per
endpoint and performs no extra material-source query or constitutive solve.

Before a proposal, the law is admitted against the immutable original thermal
source reference and accepted endpoint material. The converged endpoint source
is admitted again before acceptance. Card identity, named state, all fixed
context values and quantity types remain binding. A resolved fracture table at
one fixed temperature refuses a changed temperature, even when that change lies
inside its source validity box. Re-resolving a table at the latest temperature
does not replace original source authority. Explicit authored Johnson--Cook
coefficients may carry temperature/rate dependence across their declared domain;
this does not infer material calibration.

The returned `DuctileAccepted` contains the complete native coupled candidate
and corresponding initiation history. The caller adopts them together. Source,
temperature, history-reset, initiation, native solve and cancellation refusals
publish none of the proposed enthalpy, displacement, full plastic state, velocity
or observer state. Callers refine from the complete accepted pair to locate
initiation. Endpoint-held accumulation is first order and its crossing fraction
is plastic accumulation, not resolved event time. This gate supplies no
post-initiation degradation, fracture heat, topology change or fragments.

Five synthetic regressions in `tests/solid_plastic_enthalpy/ductile.rs` cover
independent Johnson--Cook temperature/rate utilization over successive physical
intervals, equality with the original native coupled solution and heat accounts,
initiation/refined retry/history-reset rollback, inertia and numerical-loss
separation, complete cyclic Chaboche backstress, and original-card/fixed-temperature
source refusal including attempted endpoint rebinding. All five passed against
the original native material, solid and thermal owners in the focused source
closure; this is not a full-workspace Cargo test result.

## Source-bound Timoshenko frame self-heating

`frame_vibration::thermal` and its complete exterior-fluid thermal consumer
retain native axial, symmetric bending/shear and antisymmetric bending ports.
The `fs-solid` owner supplies their common damping and original-layer work
fractions. `FrameLayerHeat` reports axial, curvature and transverse shear heat
separately; corrected shear viscosity follows its original layer area, while
curvature loss follows its original second moment. Their sum is the actual
native dissipated work, deposited once into the same finite caloric inventory.

Temperature refresh preserves original viscosity selectors, named/derived G
parents, frequency support and the authored shear correction. Geometry, density
and the actual source parents of E and G must be temperature independent in
this fixed-elasticity owner. Ordinary thermal feedback uses accepted starting
viscosities per interval; the fluid thermal owner resolves endpoint viscosity
and the complete fluid-memory trajectory jointly. Their existing work, entropy,
finite-cell, source-domain and atomic publication gates remain in force.

Live thermoelasticity retains each original Timoshenko interpolation and
consistent mass as a fixed reference Galerkin basis. The original
`ConcentricRodThermoelasticLaw` supplies affine E(T)/G(T), integrated expansion,
Helmholtz energy and entropy. Because each layer has uniform temperature and
the potential is quadratic, its RMS-curvature evaluation integrates the whole
member exactly; a separate endpoint-curvature and shear check enforces the
actual strain envelope. Named and derived G roles remain source-owned.

Each layer's symmetric spring combines bending and shear storage on the same
native coordinate. Updated viscosity, heat partition and stiffness use the
retained reference strains, while the existing complete coefficient update
checks source parameter work, Psi+TS, finite-cell heat and atomic publication.
The original complete structural basis and the same fluid/electrical/shared
assembly owners remain in force. This is a fixed-shape beam approximation:
member refinement is required when changed stiffness ratios alter the
unresolved physical fields. It does not update geometry or inertia, infer
nonlinear G derivatives, or claim exact continuum dynamics.

Focused synthetic regressions cover independent Timoshenko oscillator and
strain-work references, affine E/G thermoelastic energy/entropy exchange,
distinct layer temperatures, shear-only dissipation, full exterior-fluid
memory and warm-state refusal/retry. Native fixed-basis regressions separately
compare layer energy, loss and parameter work with integrated physical strains.


## Source-bound creep and reciprocal heat in layered shells

`shell_plastic_enthalpy::ThermoplasticShell` also accepts `CaloricCreepSource` through the shared `ThermalSectionSourceFamily`. Original standalone or joint Norton sources retain primary resistance, Arrhenius kinetics, optional sourced expansion and all native creep/plastic memories at every surface/thickness point. A standalone creep card requires no invented yield or plastic-hardening properties. Each original layer retains its own caloric cell, fixed reference mass, material identity and temperature; native conduction/observer, equilibrium, momentum and atomic feedback owners are reused.

`ShellThermoplasticEnergy::creep_heat_j` and `ShellLayerHeat::creep_j` report native creep work already included in total physical heat. Only total physical plus signed reciprocal heat enters the caloric timeline, exactly once. Numerical material/inertial loss and drilling energy remain excluded. Accepted thermal histories must match their actual enthalpy-derived temperatures, and all trials restart from the same accepted state and physical clock. Original thermal transport pinning and static-domain checks apply to creep sources as well as plasticity sources.

Two focused synthetic regressions define an independently reconstructed load/hold/unload/reload primary-creep membrane with nonzero expansion and thermal activation, late-refusal preservation, and incompatible caloric/history/card rejection. They are regression definitions, not executed results. The original fully solid, small-strain/rotation and fixed-geometry shell reduction remains bounded; measured material calibration, tertiary creep, damage, rupture and qualified service life remain separate work.

## Per-string longitudinal selection

`with_longitudinal_waves_selected` takes one `Option<LongitudinalStringResolution>`
per original string. `Some` gives that string the axial coordinates of its
resolution exactly as before; `None` leaves it transverse, with no axial
oscillators, no axial-fluid memory and no axial bearing. The calls that take one
resolution (`with_longitudinal_waves`, `..._and_loss`, `..._and_bearings`) are
the every-string-`Some` case and build the same bank bit for bit; the existing
bit-stable assertions are the regression. A selection that names no string
refuses with `longitudinal waves need a selection on at least one original
string`, and axial bearings refuse unless every original string is selected
(`axial bearings need a longitudinal selection on every original string`).

After the original string, body, cavity and transverse-gas coordinates the bank
places the selected strings' axial components in string order, then their
axial-fluid components. `StringBankLongitudinal::strings()` lists only the
selected strings; `string(i)` returns the law of original string `i` when it is
selected, and each `MaterialLongitudinalString` carries its `string_index`.
The fluid, bearing and relaxation binders address strings by that index, not by
position in the list.

The `material_strings` example reaches `None` only through the trailing word
`measured-strings-transverse` on the `longitudinal-waves` record, for
`section measured` strings. A measured whole string supplies transverse mass and
bending stiffness from its own source rows and no axial law, so it stays
transverse and is named as such; the word is off by default, and a record
without it refuses a measured string as before. The transverse-only strings'
longitudinal waves are absent from the model, not approximated.

## Per-string cylinder gas loading

`with_cylinder_gas_loading_selected` takes one flag per original string. `true` loads that string exactly as
`with_cylinder_gas_loading` does; `false` leaves it in the old sampled-resistance model: its card, modes and
damping are the unloaded ones and it has no fluid coordinates. `with_cylinder_gas_loading` is the
every-string case and builds the same bank bit for bit. A selection that names no string refuses with
`cylinder gas loading needs a selection on at least one original string`, and a selection of the wrong length
refuses. The complete mode and connection budgets count only the loaded strings' memory.

`StringCylinderLoading` gains `string_index`; `StringBankFluidLoading::strings()` lists only the loaded strings
and `string(i)` returns the loading of original string i when it has one. The fluid components follow the
bodies, one per loaded string in string order. With `with_longitudinal_waves_selected`, the axial fluid
components follow the axial ones, one per axial string that is also loaded; an axial string without gas
loading has `axial_fluid_component() == None` and no axial tail damping. The `material_strings` example
reaches the selection through the optional trailing words `strings COUNT INDEX...` of the
`cylinder-gas-loading` record; without them every string is loaded.

## Tensile envelopes with longitudinal waves and moving string ends

`SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves` carries
each source's selected tensile capacity into its original longitudinal
potential. The native `fs-nlmodal::longitudinal` strain polynomial supplies
`T_max = T0 + EA (epsilon_mean + sum_j |epsilon_j|)`. Actual oblique endpoint
extension, both transverse polarizations, their affine-chord cross terms and
the explicitly selected solid-inertia axial modes all contribute. Testing
only mean tension would miss a local tensile peak. Unretained axial channels
remain statically equilibrated as in that original reduction; this extension
does not supply their inertia or claim converged propagation.

The existing common tensile kernel combines that continuous-span axial bound
with each original layer's `E_i`, outer radius and selected static strength:
`sigma_i <= E_i (T_max/EA + r_outer_i K)`. Curvature factors are cached from
the actual transverse normalization, independently of solid axial inertia or
sparse axial mode indices. In this elastic lane the guard adds no successful-check
allocation, new force, stiffness, energy, mode, clock or solver. Native initial admission
and accepted-step gates enforce it after the original kinetic pullback;
generic copied bodies retain it and a failed trial publishes no partial state.

`SupportedStrings::tensile_assessment` diagnoses elastic candidates in the
original support coordinates without enforcing their selected strength limit.
`MaterialSupportRun::tensile_assessment` and
`ImpactSystem::support_tensile_assessment(body, wire)` read accepted native
state, including actual normalized kinetic velocities and source-bound axial
material histories when selected.
The diagnostic uses the original material-memory offsets even when other
bodies or contact coordinates precede the supported strings. Displacement-only
assessment refuses when the selected criterion needs velocity or hereditary
history. `StringTensileAssessment::tension_n` is the upper local tension bound:
elastic in the original lane, and complete equilibrium-plus-Maxwell-plus-selected
axial-dashpot tension when those source laws are selected.

`with_material_longitudinal_relaxation` retains every active arm's original
layer index, modulus increment, physical area and normalized memory address.
For each layer and retained cosine coefficient, it first sums the signed
equilibrium and actual Maxwell stresses. Only then does it take the
continuous-span harmonic envelope. Its mean includes the original maintained
preload `E_i * T0 / EA`; a branch contributes
`DeltaE_i,a * (epsilon_j - z_i,a,j / sqrt(c_i,a,j))`, where
`c_i,a,0 = DeltaE_i,a * A_i * L` and nonzero cosine channels use half that
value. Layer stresses can reinforce or cancel, and aggregate force divided
by equilibrium EA is not substituted for their load sharing.

Independently selected source axial viscosity contributes each original layer's
signed `eta_i * epsilon_dot_j` before that same harmonic envelope, with or
without Maxwell arms. Copied native bodies retain exact section, material-query
and viscosity receipts; cold binding caches each layer's viscosity. Original
normalized kinetic velocities cross the same reciprocal support-coordinate
map as displacement. The selected linear or finite-strain law supplies its
native strain rate, including the corresponding moving-end contributions.
The guarded local total-tension bound includes the same dashpot contribution.
No strain rate or layer stress is inferred from displacement increments or
aggregate damping forces.

The complete-loss source constructor binds viscosity before initial tensile
admission, so a physically relieving signed dashpot stress is retained. The
existing mechanical gate handles viscosity-only wires without creating any
relaxation state; actual Maxwell histories continue through their original
memory gate. A public empty-arm declaration remains a physical no-op and
cannot suppress initial or accepted-step tensile admission. Raw aggregate
viscosity without original layer authority still refuses for a selected guard.

The same common tensile kernel adds each layer's elastic outer-fibre bending
stress and compares the result with its selected static source limit. Guarded
axial memory is checked in the original memory acceptance gate, after actual
candidate histories are available, without an earlier elastic-only tensile
refusal. Initial and failed trial admission preserve the original transaction
semantics. The layer projection adds no physical force, stiffness, energy,
state, clock or integration method. Public aggregate-EA arms without original
layer authority remain refused for selected tensile criteria.

The envelope can be conservative when axial and bending maxima occur at
different positions. It remains a fixed-condition normal-stress estimate,
with no viscous bending stress, shear failure, fracture evolution or guarantee
between accepted temporal samples. Source bending-memory variants still
require their actual bending-history stresses. No rate-dependent strength
calibration or additional failure law is inferred from the selected viscosity.

G0/G1/G3 cases in `tests/material_string_nonlinear/tensile/supported.rs` and
its `hereditary.rs` child specify independent sparse-wave and endpoint-extension
formulas, two-plane rotation, bonded-layer load sharing, relaxed and unrelaxed
memory stresses, continuous-span containment, original body/history offsets,
source mismatch refusals, and native/prepared rollback on an actual Maxwell
tensile crossing that the elastic-only criterion would miss. The native memory
strength unit regression specifies signed coefficient cancellation before the
spatial bound. The `hereditary/viscous.rs` cases additionally cover signed
viscous cancellation and reinforcement, both selected loss laws with and without
Maxwell, initial relief, support velocity, preceding-body offsets, raw-source
refusal, and native/prepared rollback on a rate-driven tensile crossing. They
also verify that empty-arm attachment cannot remove the viscosity-only gate.


### Coupled sustained plastic/creep shell contact and finite interface heat

`shell_plastic_enthalpy::contact` advances the same source-bound layered shell, native sustained normal/Coulomb contact, optional consistent-mass momentum and original finite enthalpy network at one shared endpoint. The existing `ThermalSectionSourceFamily` preserves complete plastic, cyclic and creep histories; the common section audit and exact native resultants determine bulk physical/reversible heat, including creep only once within physical heat. The existing friction-heat network solves each massless common interface temperature and both finite-body heat transfers from original conductivity/geometry, with no prescribed heat share or added interface storage. Endpoint friction coefficients and transport refresh at the actual candidate temperatures; the callback must not hide unresolved pressure/rate/history dependence.

Every mechanical witness maps once to the original exposed layer and a distinct rigid-obstacle thermal cell. Reference skin/face positions, inward force normal, prism thickness, area and original inventory must agree; fixed footprints must remain closed along the admitted affine gap path. Pressure is native solved normal force divided by explicit area, never inferred from the shell's plane-stress normal stress. The whole-system gate closes material/caloric/kinetic storage, actual support/dead-load/obstacle work and outside heat. Material, drilling, inertial and normal-approach numerical losses remain separate from physical heat. Thermal roots stage all histories and receipts, retain the exact final mechanical response, and publish only a complete balanced candidate. One bounded cache reuses exact repeated layer enthalpies and Coulomb coefficients while source/transport receipts always refresh.

This bounded path supports fixed reference skin patches and prescribed obstacle motion; opening/changing footprints, evolving contact resistance, general impact events, wear and flash temperatures are outside it. Scheduled thermal events require interval splitting. Three authored native regressions cover independent finite-cell heat partition with temperature-dependent friction and bulk yielding, consistent-mass dynamics with numerical loss excluded from heat, and late source refusal preserving accepted states. Formatting and source/API/energy review passed; native execution remains pending the already observed shared dependency build resource blocker.

## Recovering primary resistance in solid and shell heat feedback

The existing generic `CreepSolid` and `ThermoplasticShell` owners accept
`RecoveringNortonSource`, including sourced expansion through
`ExpandingCreepSource` and `CaloricCreepSource`. Each material point retains
permanent signed creep, its nondecreasing accumulated path and a separate
recovering resistance. Temperature resolves both original Arrhenius laws over
the same physical interval. All thermal proposals restart from the complete
accepted histories; a source resolution or inspection cannot recover resistance
or advance creep by itself.

Physical creep work enters caloric storage once. Unloaded resistance recovery
alone supplies no heat; reciprocal thermoelastic heat and material/inertial
numerical losses retain their original partitions. Two manufactured native
load/hold/unload/dwell/reload regressions use independent scalar references to
check deformation, both histories, endpoint temperature, complete source
receipts and atomic refusal. The shell case includes nonzero sourced expansion.
These regression definitions have not yet run on a prescribed build runner;
their independently checked scalar trajectories are not production test results.

## Reciprocal caloric feedback through finite ductile softening

`solid_enthalpy::finite_ductile::FiniteDuctileThermalSolid` binds the original
finite ductile tetrahedral specimen to its original fully solid caloric cells.
`ReferenceCaloricMaterial` shares the existing finite thermoplastic source
admission: immutable same-card caloric object, complete typed fixed context,
original pinned Cp/density receipts, positive heat capacity and an explicit
stress-free zero-pressure caloric reference. Original reference volume and
sourced reference density must reproduce cell inventory. Optional native
consistent inertia uses that same density once; current temperature, damage
and deformation never replace the mass inventory.

The native nonpotential finite solve retains its original material, plastic
metric, initiation/damage history and velocity through every coupled proposal.
The existing finite thermoplastic Picard/IQN-ILS driver advances the original
thermal timeline, observer, events and transport. Actual caloric endpoint
temperatures and the SAME physical interval reach initiation and momentum.
Each proposal starts from the same complete accepted state. Coupling secants
belong only to the current interval and cannot advance physical history.
Source-backed initiation must match the original caloric material reference;
a calibration fixed at one temperature refuses an unsupported heated endpoint.
Authored numerical initiation retains its explicit domain without acquiring
material-card authority from the caloric source.

Compact element receipts come from the existing final constitutive audit, with
no repeated return mapping or subtraction of large cumulative dissipation.
Explicit fractions partition ALREADY-ACCOUNTED plastic dissipation and damage
release between caloric heat and two retained nonthermal irreversible energy
inventories. These are not Taylor--Quinney fractions of total plastic work;
native energetic hardening is already mechanical storage. Unheated energy is
neither discarded nor claimed to be measured fracture-surface energy, and no
later release law is inferred. Signed material endpoint defects and numerical
kinetic loss never become heat or nonthermal physical energy.

Independent plastic/damage partition checks precede the shared temperature,
sum-of-absolute-cell-heat-discrepancy, caloric-balance and combined-work gates.
The combined account retains mechanical storage, caloric storage, nonthermal
energy, kinetic change and both numerical accounts without charging physical
dissipation twice. Only the complete mechanical, thermal, observer and
nonthermal state publishes. Refusals preserve the original state for retry.
Initial binding refuses prior plastic/damage loss without its original coupled
history; subsequent intervals continue from the returned complete owner.

This extension uses the supplied fixed-coefficient finite J2 backbone, with
temperature-dependent initiation. It adds no thermal expansion, native
coefficient softening, mechanical entropy, evolving Fourier geometry, phase
change, crack topology, contact or fragment transfer. Reference-network
conductances and the heat fractions remain explicit model inputs. Coupling
convergence and first-law closure do not certify temporal or experimental
fracture accuracy. The focused definitions in
`tests/finite_plastic_enthalpy/ductile.rs` exercise independent scalar heated
initiation, subsequent release of stored energy, retained nonthermal energy,
held replay, simultaneous free-coordinate motion/heat at two physical rates,
original-source/domain/retry and cancellation refusals.

The five new coupled ductile regressions pass in the selected native execution
closure. Six existing thermoplastic regressions also pass; five older selected
cases fail identically when the original `d624b404` coupling/material/feedback
bodies are rebuilt against the same native dependencies (four native Newton
stalls and one caloric-density fixture assertion before coupling). The full
original finite ductile solid integration target passes all twelve tests.
These results do not claim a full-workspace Cargo or constellation pass: the
ordinary Cargo attempt was killed in the pinned upstream execution dependency,
and the constellation preflight reported pre-existing NumPy/pandas checkout
mismatches. No native physics, solvers or cancellation owners were substituted.

The shared source/feedback owner also preserves the concurrent finite porous
thermal adapter. Reconciled native material/solid/coupling production bodies
compile, and all five ductile thermal tests still pass. Three upstream porous
regressions pass; the fourth adiabatic-crush Newton stall reproduces with a
byte-identical residual history using the unmodified `23de8912b` upstream
feedback/material/porous bodies and identical native dependencies. No observed
new failure is attributable to this integration; the existing failures remain
explicitly outside a passing-workspace claim.

### Simultaneous recovering creep and plasticity in thermal solids

`RecoveringNortonJ2Source` composes with the existing `ExpandingCreepSource`
and native `CreepSolid` without a second material clock or projected history.
The complete fifteen-claim source and both permanent strain fields remain with
each point while primary resistance recovers independently. Physical heat is
exactly the native creep-plus-plastic dissipation; reciprocal expansion heat,
`T` times mechanical entropy storage and numerical loss stay separate. Neither
source refresh nor an unloaded recovering dwell invents heat or resets plastic
hardening. All candidate thermal roots reuse the complete accepted state.

A manufactured load/hold/unload/heated-dwell/reload regression independently
solves primary flow, plastic consistency and the thermal first law, checks all
three histories and retained receipts, and verifies proposal evaluation leaves
accepted cells, solid solution and observer unchanged. Native Rust execution
remains pending; the independent scalar reference is not a production test run.

## Source-owned body temperatures in native solid friction and heat

With `experimental-solid-contact-heat`,
`friction_heat::solid::propose_solid_contact_heat_with_body_temperatures`
maps each additional declared friction-source temperature axis to side A or B
of its original finite thermal route. Every joint trial queries the actual
endpoint caloric cell temperature while retaining the solved common interface
temperature, native pressure/speed/distance history and original source claims.
The one-result mechanical cache includes every queried interface/body
temperature; body heating or cooling cannot reuse an older friction response.

Missing, duplicate, foreign-node, primary-axis, nonthermal or undeclared
temperature ownership refuses. The original entry remains available with no
additional mappings. Fixed-footprint plane-strain mechanics, physical depth,
deposited heat, exported work and the original energy gates retain their
existing owners. No prescribed body heat split or interface storage is added.
The extended native regression covers independent two-temperature friction,
both side assignments, fixed-interface-temperature cache liveness and physical
work/heat balance. Five exact-source admission checks passed; the full native
regression awaits the pinned sibling dependencies and prescribed build runner.

## Anisotropic thermal plastic shells

The original `shell_plastic_enthalpy` family boundary accepts source-bound
thermal Hill laws with complete local-axis histories. Layer equilibrium,
reference drilling/inertia, source-resolved thermal feedback and conservative
cell heat remain owned by the existing native paths. No J2 yield/hardening
claims or state projection are required for an explicitly sourced Hill card.
Physical yield heat enters the original caloric graph once; storage changes,
mechanical entropy, reversible exchange and numerical losses stay separate.

`tests/solid_plastic_enthalpy/shell/hill.rs` exercises heated rolling and oblique
sheets under two loaded intervals and actual free unloading. An independent
uniaxial/caloric oracle checks directional flow, softened reference yield,
retained plastic memory, transverse contraction/shear, endpoint temperature,
physical heat and permanent dimensions. Source fixtures are authored numerical
cases. Scope remains fully solid, fixed material metric/frame and small strain
and rotation; the original source, geometry, step and energy refusals apply.

## Magnetized plasma absorption on finite electron/heavy gas energy

With `ionized-wave-heating`, `MagnetizedIonizedGasHeating` binds the original
magnetized half-space slab absorption to the exact original frozen gas parcel.
Its physical area/depth volume must equal the retained gas volume to roundoff;
no density, area or heat rescaling is inferred. The parcel is explicitly
uniformly mixed, while its original wave absorption is integrated exactly.

The native `fs-gas::chemistry::energy` owner receives the two signed carrier
powers once over a declared positive held-field interval. Original species
calorics supply energy differences and heat capacities; every present neutral
thermal target comes from the SAME retained electron momentum cross section.
Effective scalar electron drag alone is refused because it cannot identify
target-resolved energy transfer. The original Coulomb law supplies baseline
thermal exchange, while its signed drift correction is already in the wave
pool powers. Internal exchange cancels; coherent carrier kinetic storage,
background field storage and gyroscopic reactive power are never heat.

Native adaptive implicit exchange, source/domain and energy admission, and
physical versus numerical entropy accounts remain intact. Caller-declared
temperature increments bound the held-source step. Publication refreshes every
original carrier source and the original bound-charge card at actual Th/Te/p,
preserving immutable claims and source limits; failure or cancellation returns
no partial thermal/material result. The returned transport and background can
seed a new Maxwell solve with its own complete-field drift checks.

This is a fixed-volume, frozen-population, cycle-averaged interval with endpoint
source refresh. It does not assert continuous implicit field/thermal feedback,
source validity along an unexamined material path, a second film interface,
species evolution, spatial thermal conduction, flow or transient carrier
energy dynamics. Five authored regressions cover Maxwell-to-caloric closure
and a new field solve, signed Coulomb allocations and internal exchange,
zero-field thermal relaxation, original source/volume/domain binding,
temperature limits, cancellation and refusal of scalar-only electron drag.

`advance_reacting` uses the same held Maxwell powers with the existing
`fs-gas::chemistry::ionization_network` owner. Source-bound electron-impact
ionization, three-body recombination and spontaneous radiative recombination
share the original electron/heavy energy pools. Formation energy is stored or
released by the reaction owner, never added as a second heating source. The
coherent incident field supplies collisional work only; no separate ionizing
photon beam is inferred. Emitted recombination photons escape under the native
optically thin closure and remain a distinct energy account.

A complete `IonizedGasTransportBank` must reproduce the wave's current carrier
sources, Coulomb law and constitutive limits exactly. Its original neutral
cross sections supply thermal exchange, including any neutral target that
chemistry can create. Native reaction quadrature, adaptive time stepping,
caloric admission and conservation checks remain intact. Endpoint refresh uses
actual amounts, both temperatures and pressure, with no equilibrium reset,
population floor, inferred collision law or second charge-carrier model.

The caller bounds both pool temperature changes and every species' amount
change by `absolute_mol + relative * original_amount_mol`. Thus trace-carrier
growth is not normalized by the much larger neutral population, and initially
absent species require an explicit absolute allowance. These are endpoint
held-field limits, not error bounds on intermediate transport or proof of
continuous optical feedback. The original ion-drag and bound-charge material
cards must still apply to the resulting mixture; coordinate refresh does not
establish an unprovided composition dependence. A new Maxwell solve must admit
its own complete-field drift domain. Exactly neutral endpoints return no
free-carrier conductivity rather than a numerical floor.

Reaction evolution, increment checks, and all next-source evaluations publish
atomically. The result retains its original wave and complete transport bank,
native network result, refreshed transport/background, and separate Maxwell
work. Focused G0/G1/G3/G4 tests reconstruct monatomic sensible and ionization
storage independently, exercise the next Maxwell solve, distinguish signed
Coulomb drift work from internal exchange, permit zero-field chemical
evolution, and cover population/temperature/source/cancellation refusal. No
spatial reacting flow, line/free-free loss, diffuse recombination transport,
photoionization of the coherent wave, material-path validity or transient
field-storage claim is added.

## Two-axis string bridges in omitted-mode contact response

Each material string's original plain bridge plus one independently directed
bridge shares the full two-support omitted-mode Green matrix
`G_ab = (d_a dot d_b) G(s_a,s_b)`. Nonorthogonal and parallel axes, distinct
stations, and distinct receiving bodies are admitted. Static spring reduction
is `K_eff = (I + K G)^-1 K`; viscous reduction remains the first-order
low-frequency approximation `C_eff = T^T C T`, `T = (I + G K)^-1`. This does not
resolve omitted resonances or a full frequency-dependent viscoelastic tail.

Positive rank-one factors use at most three existing signed modal connection
columns, with the original native connection/setup budgets. There is no new
integrator. Rank-deficient compliance and zero spring/dashpot coefficients need
no inverse of G or stiffness/compliance floor. Shared double-double arithmetic
retains a weak conditional direction; only negative correlation roundoff within
128 machine epsilons is projected onto the passive boundary. The original
orthogonality admission tolerance (1e-12) and scalar/orthogonal arithmetic remain.

Fixed residual bow and stop maps use the same elastic factors in the complete
Schur complement.

Fixed residual contacts retain the affine displacement caused by the original
bridge rest extensions. The same elastic Schur factors supply
`offset_i = g_i^T K_eff r`, including scalar and orthogonal supports as well as
the full coupled two-support case. Opposite contact normals reverse this shift;
the complete reciprocal force columns and residual compliance are unchanged.

A normal stop lowers its native reference clearance to
`gap_effective = gap_original - offset_i`. This is the coordinate substitution
`U(N^T q + offset_i - G R - gap_original)`, not a new contact law. The original
material-bank stop retains its authored clearance, coefficients, quadrature
weight, physical loss and provenance. Native potential/secant evaluation,
penetration/material admission and released static preload use the same
lowered law. Source binding after preload checks the complete accepted closure,
including the omitted modes' displacement from the accepted reactions. It
cannot validate retained displacement alone or read disposable trial reactions.
Negative effective clearance represents initial interference and retains its
contact energy; it is never clamped away.

A fixed offset cancels analytically from bow/shear-film slip increments. It
adds no artificial initial shear strain, dissipation or actuator work. The
native network retains the full bridge/rest elastic energy, and omitted-tail
storage and numerical loss keep their existing accounting. Zero offset keeps
the original law and clearance bits. Rest extensions remain fixed source data;
time-varying support rest, migrating residual maps and existing incompatible
nonlinear/history lanes remain unsupported.

Resolved torsional termination springs/dashpots may coexist with scalar or
two-polarization bending residual supports, in either binding order. The actual
source-derived angular row has an exactly zero transverse prefix; rotational
lowering checks that invariant explicitly. Eliminating omitted bending modes
therefore leaves the original torque law, angular receiving-body map, rest angle
and resolved rotary inertia unchanged. The complete native modal operator still
exchanges twist/transverse energy through a shared body.

The original work-conjugate lowering remains `x = r theta`, `F = tau/r`,
`k_x = k_theta/r^2`, with the material section's actual outer radius. Native
spring/rest energy and viscous torque loss remain with the original network.
No omitted torsional compliance or extra torsional loss is inferred from the
bending Green function. This does not add a six-DOF or geometrically exact rod
termination, contact-foot separation, or dynamic unresolved poles. Bow residual
corrections leave the retained twist row unchanged; normal stop corrections
retain zero twist participation. Every torque law consumes its own native
connection column, alongside the coupled bending factorization.

More than one extra directed bridge per string, constant-loss residual bridges
and the existing mixed-contact exclusions remain outside this lane. Preparation
occurs once at binding, with no additional per-sample projection or integration
path.

Focused execution: 54 actual-source angular/torque/Schur/column cases at 4/16/64
retained bending modes, three support-axis relationships, two receiving bodies
and signed rest angles agreed with independent 100,000-mode bending tails and
the original physical angular work/energy. Full bank binding-order, torque-driven
release and bow/stop continuation regressions are authored in
`tests/material_string_torsion/polarization/rotational_residual.rs`; their native
integration execution remains unperformed in this resource-constrained session.

Support corrections retain every original receiving body. The generic friction
and normal-contact owners accept one sparse collection of signed attachments
per contact, with unique ascending original component indices; same-body rows
are combined by the material producer. Their sum enters the complete reciprocal
force/displacement column before any joint response or cross-compliance is
computed. Empty collections add no support. Total support-map work is bounded
by the original retained mode count and occurs only at binding. The actual
opposite contact surface remains independent and retains its reaction. The omitted bending
response remains quasi-static, and its numerical loss stays separate from
physical dissipation. Existing restrictions on moving contacts and incompatible
nonlinear/history mechanisms remain in force.
`tests/material_string_torsion/polarization/residual.rs` specifies a bowed X/Y
mirror with energy checks, an independent omitted-mode sine-series comparison
for oblique bow cross-compliance, active oblique-stop rotation covariance, and
equivalence when independent body coordinates are stored together or separately.
The existing material shear-film and coupled-stop residual tests independently
check each support's signed reaction, physical contact work and energy balance.

Three focused mathematical tests execute the actual two-support condensation,
including rank-one compliance at stiffness 1e30 and a weak nonparallel direction.
A bounded actual-source probe also compares 120 coupled cases against an
independent 100,000-term sine tail and matrix inverse: passive elastic/damping
operators, signed native columns, contact Schur maps and virtual work agree
within their recorded numerical tolerances. New complete-bank and contact
regressions cover nonorthogonal supports and native energy/storage behavior;
those material-admission and dynamic integration tests have not been executed.

## Oriented surface traction on two-plane material strings

`StringForceCoupling::SurfaceAlong` and
`NonlinearMaterialString::surface_shape_along` use the same source-derived
surface map. A unit force direction d and unit outward contact radius n are
expressed in the original transverse frame. The conjugate velocity is
`d dot v_centre + r*(n_x*d_y-n_y*d_x)*theta_dot`; its transpose applies the
same physical force and signed torque in the existing mechanical owner.
Positive twist rotates the first transverse axis toward the second.
The material section supplies r and rotary inertia. Measured whole-string
torsion instead retains its original measured surface-displacement basis and
supported twisting span. A radial load has exactly zero twist participation.
An eccentric load refuses missing torsion, and second-plane translation refuses
a one-plane source. No missing motion is silently projected away.

These projections are prepared once with the existing force/coordinate budgets.
The original bow-facing `Surface` behavior is unchanged. The material_strings
frontend exposes `surface-along DX DY NX NY` alongside its existing force ports.
G0/G3 regressions in the existing string-torsion and nonlinear-string batteries
cover analytical radius/inertia normalization, opposite-side torque, radial
loading, virtual work and simultaneous rotation of both transverse planes and
all physical directions. The model remains small-angle circular-section
surface kinematics; no evolving contact angle, cross-sectional deformation or
additional contact/friction law is inferred.

## Nonlinear port harmonics and sidebands through physical BEM receivers

`phs_orbit::acoustic::PortAcousticSpectrum` connects the original nonlinear
port-spectrum and infinitesimal sideband observers to the existing prepared
BEM pressure/particle-velocity banks. Ordinary harmonic peak phasors retain
their negative-time convention. Positive-time sidebands fold at their actual
signed carrier-plus-orbit frequencies: positive bands conjugate, negative bands
do not, and coincident physical lines combine coherently before source mapping.
DC remains explicit in the original port units, without an invented static
acoustic solution. No amplitude threshold removes a retained harmonic.

The caller supplies each exact-frequency source-by-port map in declared SI
source units, with its physical sign and phase. All banks must share the given
receiver identities and geometry. The returned `PreparedHarmonicField<3>`
receivers use the existing finite-window interference/intensity and calibrated
PCM owners, retaining inharmonic beats and the original source clock. Preparation
and propagation admit complete count/work/storage budgets and cancellation;
they perform no new field solve. Native BEM identities and producer uncertainty
remain borrowed in the result. This is observation of the supplied periodic
regime or small probe, not radiation feedback, an attack transient, source-unit
inference, stability proof, or a new source/spatial accuracy certificate.

`tests/phs_harmonic_acoustic.rs` exercises an original cubic-storage output
with DC/even/odd harmonics through actual incident-plus-scattered BEM fields
and the PCM stream, conjugate sideband collisions, independently integrated
finite-window beats and vector intensity, static-only output, and physical
binding/budget/cancellation refusals. Forced and autonomous cycle-work audits
also refuse grids that cannot resolve the authored input or retained state
harmonics; nonlinear power still needs successive-grid convergence evidence.


## Exact compact radiation from linear string harmonics

The fixed-end linear string renderer and its material, thermal, concentric and
finite-gas runtimes evaluate compact-strip pressure from the accepted modal
state and authoritative damping rates. Under held force and coefficients,
`a = f - omega^2 q - d v` and `j = -omega^2 v - d a`. The original signed area
and first-moment weights retain odd-mode acceleration and even-mode jerk
radiation. The observer adds no allocation or sampled derivative history, so
even partials avoid backward differences' attenuation and half-sample phase lag.
Material and gas updates evaluate the same expression under their final held law.

Constrained moving-end interval-mean observations retain their clock. This
remains a compact strip observer without radiation loading or common
propagation delay; ideal force, release and parameter-jump impulses require a
resolved transition and are not
reconstructed. Analytic high-band amplitude/phase, two-grid public-runtime,
damped held-load and material/gas-update regressions cover the changed behavior.

## Complete nonlinear string endpoint radiation

Fixed-end Kirchhoff--Carrier strings now evaluate acceleration and jerk at the
proposed accepted endpoint using the original full pHS law. With constant
structure and held input, `xdot = (J-R) grad H + G u` and
`xddot = (J-R) D(grad H) xdot`; canonical mass-normalized momentum rows supply
the observer's modal acceleration and jerk. The storage-owned continuous
tangent retains shared cubic stretching, bending, elastic distributed contact
and all Maxwell memory coordinates. Four reusable state buffers and one zero
input buffer are prepared once; this observation and these admitted borrowed
storage actions allocate nothing while running. Pressure admission precedes
mechanical publication and preserves the existing texture rollback.

This differentiates the physical storage at the numerical endpoint, not the
numerical step map, and does not make nonlinear propagation exact. Existing
bow, tangential-friction and Hunt--Crossley input forces are held over each
step: their continuous feedback derivatives and ideal jump impulses are not
claimed. Alpha=1 contact at exactly zero penetration and nonzero normal
velocity has no classical two-sided effort derivative and refuses explicitly;
quiet contact, active/separated contact and alpha>1 entry retain their actual
derivatives. No finite-difference fallback changes that meaning. Joined
moving-end reactions remain interval-mean quantities. Compact-strip geometry,
sampling/aliasing limits and the existing polarization model are unchanged.

Independent force/jerk regressions cover bending dispersion, shared stretching,
modal damping, held-load changes, mixed active/inactive elastic contact, and
radiation caused by memory relaxation with zero instantaneous modal velocity
and acceleration. A renderer regression checks amplitude and phase of an even
partial at 80% of Nyquist; contact-corner coverage checks the named refusal.
These native regressions require execution through the real dependency graph.

### Coupled rod fatigue with bonded layered transverse shear

`RodMotionFatiguePath::new_with_layered_transverse_shear` selects the original
solid observation owner's bonded common-Poisson circular recovery for native
motion. Distinct original layer moduli and source identities reach the same
S-N and Findley tensor histories. The wrapper retains the existing mechanical
trajectory, accepted-pose observer, work and rollback owners; it supplies no
new integrator, averaged material law or additional mechanical energy.

### Native sustained finite-solid friction heat with source feedback

`friction_heat::finite_solid::propose_finite_solid_contact_heat` closes endpoint
body/common/declared-mean source temperatures, the original sourced sustained
sphere/finite-solid trajectory and the existing finite enthalpy network together.
Each original feature retains its source binding, ordered materials, caloric
objects and bounded reference thermal patch. Only unique loaded, nonopening,
passive midpoint tangential impulse work is eligible for declared heat
conversion. Exported tangential work, signed normal/contact remainder and native
integration defect stay separate in the combined first-law balance; no residual
or normal loss is deposited as friction heat.

Every trial restarts the accepted solid, rigid bodies and original source history.
Exact repeated full source queries reuse one bounded mechanical result while
transport refreshes; common interface temperature alone is never a cache key.
Source, geometry, authority, thermal, clock or balance refusal returns no partial
accepted state. Fixed witness neighborhoods and ordered face normals are checked
at native source evaluations and the final endpoint. This does not add a swept
thermal-footprint certificate. The joint candidate retains original synchronous
endpoint times and all mechanical/history/thermal reports.

This is fixed-reference-patch interface temperature feedback with frozen bulk
hyperelasticity. Pressure-dependent source coordinates, moving footprints, bulk
expansion and impact heating remain outside the model. The existing trajectory
and thermal interval still need time refinement; nonlinear convergence is not
time-accuracy evidence. The native regression definition in
`tests/finite_solid/contact/step/sustained/source/heat.rs` compares an independent
body-temperature-dependent friction/impulse/caloric solution, zero conversion,
source path, fixed-patch refusal and late source rollback. Isolated execution of
the production work-account function with metadata fixtures validates its heat,
export and remainder partition, not the complete native multiphysics solve.

## Cylindrical and spherical moisture diffusion

`moisture::radial::MoistureRadial` advances bound-water diffusion on nonuniform
cylindrical or spherical shells using their actual dry-skeleton storage volumes
and exact radial face resistances. Original source-integrated Kirchhoff laws
retain moisture-dependent and degenerate dry-front transport; finite sourced
surface exchange uses each boundary's curved area. Conservative backward Euler
uses a reusable native tridiagonal workspace, bounded nonlinear iteration and
independent endpoint residual/water-balance checks before publication. A solid
center is sealed. This is a fixed, homogeneous, isothermal radial reduction;
axial/angular gradients, arbitrary tensor orientation and swelling feedback
require their respective spatial/coupling owners.
Analytic cylindrical/spherical steady potentials, a spherical transient
eigenmode, nonlinear dry-front balance and finite surface flux define the focused
verification in `tests/moisture/radial.rs`; they are not measured validation or
a continuum-error certificate for an arbitrary mesh or time step.

### Prescribed temperature histories in radial diffusion

`MoistureRadial::at_prescribed_temperature` re-evaluates the original pinned
multiaxis diffusivity at every occupied shell content while retaining each
shell's dry mass, water, geometry, accepted time and original density/reference
identity. A supported enclosing temperature/moisture box cannot override the
actual source triangulation. Degenerate dry fronts retain the original
Kirchhoff path, and an unsupported cell or cancellation refuses the whole event.

Radial backward Euler, SDIRK2 and finite-atmosphere admission consume the current
transport temperature. Source-qualified sorption and finite surface laws must
match that temperature; the cold density datum does not silently drive a warmed
boundary. Finite solutions already consume the same current inventory context.
Transport intervals must end at externally prescribed temperature events.
This path supplies no heat balance, thermal expansion or temperature evolution;
the native mechanical owner still refuses an unqualified thermal strain.

### Source-partitioned coatings and curved material layers

`moisture::radial::layered::{MoistureRadialLayer, MoistureLayeredRadial}` joins
contiguous native cylindrical or spherical shells while retaining each original
dry inventory, nonlinear diffusivity and same-card monotone sorption source.
The shared layered activity/storage/Picard kernel serves both slabs and curved
stacks; each geometry supplies its actual face paths. Radial layers use exact
cylindrical/spherical resistance integrals, unequal shell masses and the actual
curved interface and exterior areas. Existing slab deformation-metric feedback
continues to resolve the proposed state at each nonlinear assembly.

Ideal contact equalizes activity while allowing material-specific dry-basis
moisture jumps. Explicit finite activity conductance permits an interface
activity jump with one equal/opposite water transfer. A zero-area center stays
sealed. Source/context, geometry, temperature, clock and support mismatches
refuse; native source/range/work/residual/water checks precede publication of
the entire stack. Geometry and temperature are held fixed; this rung does not
add hysteresis, moving interfaces, sorption heat or liquid capillary flow.
Five authored physical regressions define curved steady partition/transfer,
unequal-density closed transients, homogeneous-owner reduction, nonlinear dry
core/coating uptake and transactional geometry/source/admission refusals.

### Finite solution activity on radial geometry

`FiniteMoistureSolution::advance_radial` and `advance_shared_radial` connect
native cylindrical and spherical diffusion to the original finite aqueous
inventory. Each shell contributes its actual dry-mass-weighted water change;
one common activity root retains all nonvolatile solutes and any admitted
neutral/ionic or kinetic phase owner. Curved half-cell resistance, nonlinear
Kirchhoff transport and sourced surface exchange remain in the radial solver.
Stable specimen IDs fix signed transfer summation, and all fields and the
solution publish together after the original water and source checks.

The corresponding `advance_radial_sdirk2` and `advance_shared_radial_sdirk2`
solve activity at both physical stages. Both start from the original solution
inventory; the final report already contains the full weighted water transfer.
Kinetic crystallization explicitly refuses this second-order path because its
endpoint-frozen phase law supplies no RK stage equations. Equilibrium crystal
stocks remain finite and may buffer the activity. Smooth-branch temporal
accuracy requires refinement; phase events need their own resolution.

The coupling retains prescribed common temperature, equal pressure and an
equilibrated vapor path with negligible vapor storage. It does not reinterpret
vapor transfer coefficients as liquid-immersion permeability or move solutes
into the solid. Five focused authored regressions define analytic curved uptake
and release, concentration feedback, unequal-shell shared conservation/order,
second-order convergence, finite crystal buffering and transactional refusal.

### Second-order heterogeneous water storage and interface transfer

`MoistureLayeredSlab::advance_sdirk2` and
`MoistureLayeredRadial::advance_sdirk2` use Alexander SDIRK2 through the same
native layered storage/activity kernel. Both stages retain each original
nonlinear sorption and Kirchhoff diffusivity, actual dry-cell inventory and
explicit interface/exterior law. Stage two starts from the original water
storage and receives the actual first-stage face rates as physical water loads;
neither an activity increment nor a frozen capacity substitutes for stored water.

`MoistureLayeredSdirk2Report` separates the first diagnostic stage from the
complete endpoint. Exterior and material-interface water transfers in the
endpoint already include both RK weights; surface activities and source
receipts describe the final stage. Shared source preparation, native
tridiagonal iteration and admission remain in use for both spatial geometries.
The existing backward-Euler and deforming-slab paths retain their original
behavior. This SDIRK entry point holds geometry, temperature and reservoir
activities fixed. A nonmonotone stage returns typed `ActivityStageRange`,
which the existing duration controller may retry at a shorter interval.
No clipping, first-order replacement or finite-reservoir SDIRK claim is implied.

### Finite vapor exchange with coated radial bodies

`moisture::radial::layered::atmosphere` couples single or shared coated spheres
and cylinders to the original finite-vapor inventory. Each humidity trial
advances the complete original body with every material-specific sorption,
nonlinear diffusion and interface law intact. Compensated local moisture
increments, weighted by the actual unequal dry shell masses, supply one gas
water account. Canonical specimen IDs give deterministic shared summation and
result order. An interior material interface is never an atmospheric port.

Exterior source receipts resolve the pinned end-material sorption claim at the
accepted solid-side surface activity; explicit film resistance can separate it
from gas humidity. Authored exterior seals remain sealed while internal
redistribution continues. Source/context, common clock/temperature, zero-area
center and total cell budget admission precede the original water/carrier/total
mass checks and atomic publication. Temperature, free gas volume and reference
geometry stay fixed; this path supplies no condensation, carrier permeation,
hysteresis or thermal feedback.

### Finite atmosphere exchange on radial geometry

`FiniteMoistureAtmosphere::advance_radial`, `advance_shared_radial` and
`advance_paired_radial` couple the original radial fluxes to the existing finite
gas inventory owners. The actual implicit endpoint RH drives sourced sorption
and surface transfer; every probe restarts the original specimen. Compensated
local water increments use each shell's actual dry mass. Shared chambers debit
one gas inventory for all stable-ID-ordered specimens; paired chambers retain
separate signed inner/outer transfers even when total solid water is unchanged.
The original atmosphere owners resolve changed composition and pressure and
check carrier, water and total mass before all endpoints publish together.

Native solid adapters resolve compatible deformation and wet vibration inertia
after the accepted single or paired atmosphere candidate. A late constitutive
or mechanical refusal preserves every original gas and solid input. Gas
pressure is not automatically imposed as a mechanical traction. Temperature and
gas volumes are prescribed; condensation, carrier permeation and sorption heat
remain separate physics. Source-supported RH brackets and explicit work budgets
remain caller inputs. Seven radial gas regressions cover analytic exchange,
finite surface resistance, opposing chamber flow, unequal shared inventories,
dry-front wetting, sealed redistribution and atomic refusals; a native solid
regression adds finite-vapor shape/inertia and late mechanical refusal.

### Second-order radial transients

`MoistureRadial::advance_sdirk2` uses the same two-stage L-stable SDIRK method
as the existing slab and Cartesian owners. Both stages reuse the original radial
Kirchhoff law, unequal shell storage, curved surface exchange and bounded native
tridiagonal nonlinear solve. The first-stage explicit load is assembled from
actual conservative face rates, never rounded water differences. The endpoint
water account uses the RK-weighted sum of both stage boundary transfers.
`MoistureRadialSolid::advance_sdirk2` resolves compatible native mechanics after
the complete admitted endpoint, preserving atomic transport/mechanics publication.

L-stability does not guarantee monotonicity: a nonphysical stage returns the
existing typed `StageRange` refusal and requires a shorter time interval. It is
not clipped or replaced by backward Euler. Conditions are constant within the
call; callers split intervals at environmental events. Independent unequal-shell
linear modes, a nonlinear dry-front ODE, finite surface resistance and weighted
water balances define six focused native regressions. A converged stage is not
a global temporal-error estimate.

### Second-order and adaptive finite-vapor radial transport

`advance_radial_sdirk2` and `advance_shared_radial_sdirk2` on the original
finite-atmosphere owner solve the actual gas humidity at each implicit stage.
First-stage face rates supply the explicit RK load; the second stage balances
the original gas inventory and original solid storage against the full weighted
water transfer. Diagnostic intermediate gas/solid states are never debited twice.
The prepared Kirchhoff source interval retains the first stage's admitted range
for later humidity probes, while the algebraic reference shift and actual
active-face physical range checks remain unchanged.

`advance_radial_adaptive_to` and `advance_radial_sdirk2_adaptive_to` control
a complete finite-air duration using the existing temporal driver. Acceptance
checks every shell's moisture, noncancelling local water error, gas vapor mass
and RH. Accepted fine-step transfers alone enter the fixed whole-window water,
carrier and total-mass account. Source/bracket limits remain explicit; a failed
duration publishes neither inventory. Matching native-solid methods resolve the
compatible final shape and wet vibration pencil after these gates, with late
mechanical refusal preserving all inputs. Second-order endpoints also
support multiple specimens sharing a chamber. Adaptive finite-air durations
currently couple one specimen and one chamber; paired-chamber exchange retains
its backward-Euler path. Eight native regressions define
independent coupled exponential limits, nonlinear curved-film fluxes, stage
source support, unequal-shell/gas accuracy and whole-call rollback. Their error
estimates are local, not global continuum or experimental-validation claims.

### Accuracy-controlled radial durations

`MoistureRadial::advance_adaptive_to` and `advance_sdirk2_adaptive_to` reuse
the existing whole-duration temporal controller around the original backward
Euler and SDIRK2 endpoints. The shared borrowed inventory now explicitly owns
either uniform or per-cell dry masses. Generic water, caloric, reservoir and
solid accounts consume that physical measure, so nonuniform radial shells are
never interpreted as equal storage cells.

Each comparison resolves both componentwise moisture and the L1 sum of actual
shell water discrepancies; opposite internal errors cannot cancel in a closed
specimen. Relative scales use transported changes. For radial calls,
`reservoir_water_absolute_tolerance_kg` also controls this shell-water norm.
Only accepted fine-step transfers enter one fixed whole-duration water budget.
The existing controller charges rejected/stage work, bounds retries, handles
typed physical stage refusals and publishes no successful prefix after a later
failure. Matching native-solid convenience methods resolve the mechanical
endpoint after the complete admitted duration. The estimate remains local and
asymptotic; boundaries are held within the duration and environmental events
require explicit cuts. Six focused native regressions use unequal-shell modes,
nonlinear dry-front solutions, an independent mass-weighted error ratio and
work/cancellation rollback.

### Adaptive radial diffusion with a finite solution

`FiniteMoistureSolution::advance_radial_adaptive_to` and
`advance_radial_sdirk2_adaptive_to` use the original bounded duration controller
with complete native finite-solution backward-Euler or SDIRK2 endpoints.
Every trial solves actual solvent-volume/activity feedback with the original
curved transport paths. Error estimates compare each shell content, the
dry-mass-weighted L1 shell water discrepancy, finite solvent mass and solution
activity. Opposing shell errors cannot cancel, and pure-solvent uptake remains
controlled even when its activity stays one.

Only accepted fine endpoint exchanges enter the cumulative boundary account.
The original-to-final shell/face balance and combined solvent/shell water balance
must each meet one fixed whole-duration tolerance. Both SDIRK stages count
against the work budget, and any late failure discards the entire candidate.
Adaptive radial solid wrappers resolve shape, stress and wet inertia once from
the accepted final moisture, then publish the complete dependent result together.

The existing prescribed-temperature, equilibrated vapor bridge remains in use;
these methods do not add liquid immersion or osmotic traction. Native kinetic
phase evolution remains first order under BE and is refused by SDIRK2. Individual
phase-component errors and phase heat are not separately controlled. Local
step-doubling estimates do not promise a global accuracy bound.

### Radial sorption energy and finite-solution mechanics

`FiniteMoistureAtmosphere::advance_radial_with_sorption_energy` and its SDIRK2
counterpart integrate the original differential isosteric heat over each actual
shell's dry-mass-weighted loading change. The existing caloric owner closes
original gas plus bound-water internal energy against the heat supplied by the
external isothermal thermostat. Loading-dependent internal redistribution can
require heat even when no water leaves the specimen. The intermediate SDIRK
stage is never added as a separate energy debit. Source work and energy-balance
failure withhold the complete transport/energy candidate. This is global heat
duty, not a reconstructed local heat equation or temperature feedback.

`MoistureRadialSolid::advance_finite_solution` and its SDIRK2 counterpart resolve
compatible native deformation, stress and the wet prestressed vibration pencil
from the exact finite-solution diffusion endpoint. Late mechanical/source
refusal withholds both water inventories and their dependent solid state.
Osmotic pressure remains diagnostic; no immersion flow, inferred traction or
elastic-work feedback into chemical potential is added. Five authored tests
define unequal-shell heat integrals, sealed redistribution, late caloric
refusal, native affine swelling/wet inertia and late mechanical refusal.

### Radial diffusion driving native solid mechanics

`moisture::solid::MoistureRadialSolid` maps every original cylindrical or
annular water shell into a conforming three-dimensional tetrahedral mesh. Local
transported moisture selects the original elastic law and integrated directional
swelling source; the existing native solid owner resolves compatible displacement,
stress, reactions, deformed geometry and the signed prestressed K_e+K_g/M pencil.
A nonuniform water profile therefore produces compatible stress and local wet
inertia rather than independent shell expansions. Material axes are fixed in the
xyz mesh frame; radial/tangential grain rotation is not inferred.

Each shell's exact curved dry inventory is distributed over its actual polygonal
tetrahedral volumes. The user explicitly bounds the curved/polygonal volume
discrepancy; the corresponding inertia-density correction vanishes under angular
refinement. Axial, circumferential and radial mechanics still require ordinary
mesh convergence. The water metric remains on the original reference chart:
small deformation does not imply stress-assisted or spatially swelling diffusion.
This owner admits cylinders and annuli; spheres require a separate native mesh
owner. Source-domain, geometric, work, cancellation and late mechanical refusals
leave the original moisture state untouched. Geometry tests independently check
shell volumes, seam conformity and refinement; native integration tests define
affine free swelling, gradient-induced stress, exact conserved wet inertia and a
one-DOF generalized modal response using authored verification laws.

## Coating strain about a supplied static modal offset

`MaterialShellImage::prepare_coating_strain_with_offset` and the corresponding
assembly-image method admit `PreparedShellCoatingStrain::evaluate_with_static_offset`.
The supplied real modal offset and coherent complex modal amplitude use the same
original source-bound mode basis, facet geometry and coating sections. Each
reported film bound covers the full physical layer, complete triangle and
complete harmonic cycle about that offset; independently acceptable offset and
oscillation are not assumed acceptable when combined.

Preparation admits the extra real field, modal accumulation, storage and bounded
static-plus-harmonic maximization work before observation. A successful query
retains complete film coverage; invalid input, unsupported preparation and
cancellation do not consume the successful-query allowance. The reported
sufficient amplitude scale multiplies both supplied fields together. It is not
an oscillation-only allowance at fixed preload. The observer does not infer
thermal/cure eigenstrain or solve nonlinear prestress equilibrium.

Three integration regressions cover combined-limit crossing, phase invariance,
native snapshot containment and transactional admission/refusal. Their complete
source-material/modal fixture remains unexecuted in the available dependency
environment. The underlying native shell bound passed five focused numerical
tests, including complete-domain containment, degeneracies and extreme scaling.

## Original thermal hardening surfaces on reciprocal caloric solids

`FiniteThermoplasticMaterial::try_new_with_tabulated_hardening` binds an
explicit tabulated finite thermal model and its original strain/temperature
observations to the exact same-card caloric chart. The complete original
hardening source and every selected observation remain inspectable through the
resolved family. Original density at the stress-free reference temperature still
sets reference mass; expansion and density metadata cannot change that inventory.
Both scalar and tabulated constructors use the same caloric source/context and
positive-total-heat-capacity admission.

`FiniteThermoplasticSolid::try_new_with_tabulated_hardening` prepares these
families through the same cell mapping and native tetrahedral constructor.
Explicit claim selection and a per-source preparation budget are supplied by
the caller. Existing equilibrium, retained inertia, accepted-state transactions,
physical interval and reciprocal heat iteration remain the only step owners.
The prepared material can also enter the existing native finite thermal sheet
and `FiniteThermoplasticSheetCoupling`; no second heat solver is introduced.

Heating can cross a measured hardening knot and advance plastic strain before
melting. Its hardening-entropy change is retained in signed reversible caloric
exchange, independently from positive physical plastic heat. Numerical endpoint
work defects remain separate. The same accepted plastic metric, accumulations,
temperature and caloric inventory start every feedback proposal.

Two coupled regressions author an independent coaxial table/entropy oracle,
solve a free displacement with temperature, cross a hardening knot on heating,
and cool/unload to retained permanent shape. They also check original receipts,
foreign caloric rejection, bounded temperature/strain support, complete
nonvirgin state retention after refusal, and successful retry. Their independent
scalar energy identities and Rust syntax were checked; native regression
execution remains unverified in the resource-constrained session.

## Spatial pressure in moving partial-slip friction and heat

`friction_heat::moving::partial_slip::MovingPartialSlip::new_with_pressure_profile`
accepts `MovingPressureProfile` from the original normal-patch authority. Its
positive, piecewise-constant pressure cells exactly cover the unchanged strip
footprint and integrate to the original normal load, to arithmetic precision.
The original patch identity, source, authority, load and geometry must match the
retained checkpoint. Gaps, overlap, foreign receipts and load mismatch refuse;
there is no extrapolation, load renormalization or inferred pressure shape.

Pressure stays in footprint coordinates while the two native Jenkins-region
histories travel with surface A material. Exact intersections split parcels at
thermal/material and pressure boundaries without averaging strains or resetting
kinetic branches. Sorted interval cursors keep this splitting proportional to
the actual intersections. Each source callback receives its local pressure,
actual regional normal force and original caloric cell. The original return
map still owns local force, microslip, stored energy and physical heat; its
reference-normalized response is multiplied by the actual regional area share.

Stiffness, region decomposition and the storage-source identity remain fixed.
A lower local pressure may relax retained shear energy through the original
return map, even at zero additional imposed slip. Only its physical plastic
loss heats the finite cells. Numerical return-map loss and elastic energy
leaving the footprint remain separate accounts; no parameter work or normal
work is invented or thermalized. Raising pressure preserves loaded history.

Six new physical regressions cover uniform-profile equivalence, early edge
microslip with a sticking core and spatial heat, transport across pressure
boundaries with independent release/outflow balances, local pressure plus
endpoint-temperature source feedback, exact source/load/coverage admission,
and bounded transactional refusal/retry. They accompany the existing moving
brush tests for work, frame covariance, history and temporal refinement.
Pressure-space, thermal-space and time refinement remain independent. The
existing pure-shear, fixed-load, fixed-width and first-order transport/contact
composition is retained; zero-pressure gaps require a separate normal/contact
topology, and no resolved continuum traction or measured calibration is inferred.

### Evolving normal load and pressure without resetting moving shear

`MovingPartialSlip::step_candidate_with_pressure_profile` accepts the original
normal owner's next positive load and piecewise-constant pressure field on the
same footprint. The native partial-slip checkpoint admission retains exact
patch/model/source/authority and ordered-interface lineage. Semiaxes, pressure
moment, shear stiffness and its source stay fixed. New pressure intervals must
cover the unchanged physical footprint and integrate to the supplied new load;
no coefficient or pressure rescaling repairs an inconsistent receipt.

Transported material strains and kinetic branches enter the original regional
return map unchanged. Unloading can release stored shear at zero imposed slip;
only its native plastic dissipation heats the original finite thermal bodies.
Numerical return loss and elastic outflow remain separate. Increasing load
retains the surviving strains without inventing preload. Local source queries
receive the new regional pressures together with both actual caloric cells and
their solved interface temperatures. Pressure, all material histories and heat
publish atomically. `normal_patch()` exposes the accepted normal receipt, while
`reference()` retains the original construction checkpoint; ordinary subsequent
steps continue with the accepted pressure automatically.

This is a first-order endpoint-pressure split. Normal-path refinement is
independent of pressure-cell and thermal-mesh refinement. Normal indentation,
normal storage and actuator work belong to the supplied normal owner; this
brush neither infers nor thermalizes them. Separation, footprint growth,
pressure-dependent shear stiffness and a coupled normal equilibrium solve
remain outside this path.

Five G1/G3/G4 regressions verify analytic unloading/reloading, local pressure
redistribution at constant total load, convergence of physical heat under
16/32/64 load subdivisions, lineage/geometry/load refusal before source work,
and atomic late refusal, cancellation, budget exhaustion and retry. Together
with the preceding moving-memory, pressure and two-temperature regressions,
all 22 focused native tests pass. No measured calibration is inferred.

## CST/DKT hereditary material memory on the original modal basis

`plate_viscosity::shell::relaxation::ModalShellRelaxation` projects per-facet
passive excess A/B/D stiffness spectra through the original CST/DKT strain
owner. Its modes must be orthonormal in the original sparse mass, including
consistent director/translation terms; the instantaneous-viscosity mass
check is shared unchanged. The equilibrium pencil and retained coordinates
remain fixed. Source applicability, matching mesh/basis identity and the
interpretation of supplied tensors as excess above relaxed stiffness remain
explicit caller obligations.

Exactly equal relaxation times aggregate before the shared compensated PSD
factor sees each original dimensional modal Gram matrix. Different times,
resolved weak directions and all modal cross terms remain present. Each
resulting row c [1/s] carries the existing Maxwell storage `(c.q-z)^2/2`
and irreversible power `(c.q-z)^2/tau`; numerical drilling acquires no
constitutive state. Spatial work, complete factorization, branch embedding
and material-state counts have explicit checked budgets. No dense projection
or factorization occurs during playback.

`branches` embeds these rows into existing fs-phs interleaved mechanical
coordinates. `attach_to` composes them through the original
`CoupledModalSystem::with_additional_material_relaxation`: relaxed history
sets new z=c.q, unrelaxed history sets new z=0, and earlier material histories
remain unchanged. The target component must have zero intrinsic modal
damping; separately selected viscous, support, glue and radiation channels
can coexist. Existing cold-construction, instantaneous-frequency, energy and
augmented-state admissions remain binding. The existing exact propagator and
independent material-loss integral are the only time and heat owners.

Three focused native tests pass: full spatial/modal/channel storage for both
mass models and every retained time, branch embedding/refusal/zero cases,
and actual fs-phs exact propagation against independent Maxwell equations
with separately integrated heat and input work. Both initial histories pass
the 16/32/64-subdivision RK4 convergence checks and energy closure without
tolerance changes. The additional full coupled-attachment regression covers
component offsets, prior histories, initialization, loss counting and cold
binding; it is authored but unexecuted in the available dependency closure.

