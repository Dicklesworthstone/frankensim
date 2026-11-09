<!-- fs-couple contract, part 2 of 8; a higher number is newer. Index: ../CONTRACT.md. Keep this file under 900 KB: the GitHub contents API returns empty content above 1 MiB, and a writer that sees an empty file overwrites it. -->
## Source-resolved Hertz islands on shared finite thermal cells — v0

`friction_heat::hertz::multiple::propose_rigid_hertz_island_heat` couples
multiple original Hertz patches to one original native body island and one
finite enthalpy network. Each stable patch identity binds its pinned original
coefficient source, unchanged normal load/radius, one distinct original thermal
route and its own solved common-interface temperature. Explicit pressure and
the T/p-only source-axis refusals remain intact. Conductivity inputs retain
canonical original link identities; a precomputed permutation maps thermal
link order to native patch order without repeated trial searches.

One native island is prepared over the actual representable interval
`(accepted_time + requested_dt) - accepted_time`. Every nonlinear trial
refreshes all original source claims at their respective interface temperatures,
solves all mechanical patches together from unchanged accepted momenta, and
routes their actual dissipation through `FrictionHeatNetworkStep`. That original
network sums powers into shared cells once, suppresses the replaced direct
contact links once, and retains all ordinary outside heat and other conduction.
No fixed heat fraction, duplicate capacity, prescribed face temperature or
numerical-loss heating is introduced. Geometry, circular area, centre and
normal must match the original native contacts.

Only a complete body/source/thermal candidate leaves the owner. Per-patch
mechanical and interface audits remain in force, and both direct
`Delta K + Delta H - Q_outside` and the sum of the original mechanical/network
residuals must pass independently. Final-candidate caloric entropy can be
requested from the same original network. Source failures, exhausted budgets,
late audits and cancellation preserve all accepted owners and exact retry.
Allocation, unknown count, nonlinear iterations and total trials are explicit
bounds; native controls are admitted before trial allocations.

Frozen geometry/normal data, the native response-rank and endpoint-event limits,
original orthogonal finite cells, common massless faces and backward Euler
remain the approximation scope. This does not consume retained elastic
partial-slip memory or rolling loss and does not establish root uniqueness or
experimental validity. Seven focused G1/G3/G4 groups in
`tests/articulated_friction_heat/hertz_multiple.rs` use an independent
three-body/two-patch mechanical response and eliminated finite-cell thermal
root, check source receipts and entropy, canonical permutations and signed
outside heat, stationary two-route conduction, source/geometry/ownership
refusals, loaded late failure/cancellation retry and large-clock synchronization.

### Complete vibrating history at moving-end Dirac attachments

The moving-end renderer consumes the original `PlateBankImage` through
`dirac_body::join_complete_body`. Mixed compact/Vk coordinates, full material
damping, Maxwell memory and attached acoustic histories remain in their
original bank order after the complete string state; any explicit fluid load
follows them, and the existing join multiplier remains last. The body force
port participates in the original common-velocity `join_port`; its independent
signed outward-volume port couples through the existing ratio-one transformer.
The surviving fluid source port is mapped explicitly after the original string
force ports. An unloaded body's pressure port stays open with zero applied
pressure and is not mistaken for a volume-flow source. A reed or volume source
without an independent duct inlet refuses instead of being silently ignored.

Joining accepted histories requires compatible initial attachment velocities.
Admission uses the original efforts and force columns with a relative
roundoff scale from their absolute products. Incompatible moving histories
refuse; no velocity projection, hidden impact, state reset or material-memory
reconstruction occurs. The multiplier has no storage and starts as a zero
representative. Exactly declared quadratic joins prepare the shared
`PreparedCommonFlowZoh` once, retain every physical coordinate, and return the
instantaneous endpoint reaction. Original material heat and source work are
integrated independently by that owner. Its complete-state and dense limits
remain; the renderer reuses candidate-state and mean-output buffers.

Nonlinear storage keeps the original descriptor/Gonzalez update and the
original complete pole-rate bound and midpoint admission, including retained
acoustic poles. Exact held-input evolution does not make state-dependent
exciters monolithic or supply a source-bandwidth/anti-aliasing guarantee.
The existing held-input, contact and observer conventions remain explicit.

The original bank observer reads its complete physical prefix at the declared
string offset, so compact and Vk modes cannot overlap or double count in the
pressure calculation. The string contribution uses the same compact strip
observer as fixed-end strings, with the signed area and midpoint first moment
of its original half-integer cosine shapes. Both its monopole acceleration and
dipole jerk therefore survive body/duct composition. Acceleration retains the
original interval-mean momentum increment; jerk differences consecutive accepted
interval means. The first interval primes its own history without inventing a
release impulse or an endpoint reaction from the descriptor's mean multiplier.
Body/duct observation and material histories keep their original owners. The
existing positive-string-axis compact approximation and absence of string
radiation back-reaction remain. G1 checks compare sine/cosine moments with
independent spatial quadrature and actual bare/joined pressure with a damped
oscillator derived by eliminating the common-velocity constraint.

Focused native tests cover independently written
mechanical/fluid DAE evolution with 20/40/80-step refinement, work/heat and
constraint closure, nonlinear histories, malformed/incompatible joins, and
128 exact steps of accepted warm radiation memory at 48 kHz. Actual renderer
checks cover mixed and Vk-only bodies, signed areas, declared radiation at
48 kHz, render-length prefix identity, an exact quadratic 10 kHz common mode,
and retention of the nonlinear clock refusal. The native target uses the
actual exercised physics; unavailable bow/texture and cavity-source setup
are explicitly outside its tested scope. Full workspace and measured
real-time qualification are not asserted.
A dependent two-body/two-patch fixture also couples a shared native support
to two distinct original finite thermal pads carrying that same body's identity;
its independently derived full mechanical response and finite-cell root check
both body momenta, each physical heat contribution and all temperatures.

## Native circular wetting between unlike surfaces

`contact::capillary::AsymmetricCircularContactLineCapillaryModalSystem` binds
the original axisymmetric bridge with two independently evolving circular
footprints to the shared native mechanical transaction. Both ordered source
identities, surface histories, wetting thresholds and line resistances remain
with the original liquid owner. The actual body gap drives the joint footprint
solve; one line may pin while the other moves. Geometry and line perimeter are
three-dimensional, and each line exposes its own viscous and pinning losses.

Mechanical probes all start from the same accepted liquid state. A reserved
final replay evaluates both footprints on actual represented body motion,
with original component, connection, external-force and prescribed-gap work
checks. Geometric and quadrature uncertainty must fit the native compatibility,
force and energy gates before bodies, liquid history, diagnostics and clock
publish together. Refusal and cancellation leave both original footprint
histories and the native states retryable.

The owner inherits the original connected, quasistatic, acute-angle,
nonpositive-pressure, interior-neck and low-speed applicability limits. It
does not add bulk liquid inertia, evaporation, topology change or a full 3-D
stability claim. G1/G3/G4 native tests cover unlike-source pin/move dynamics,
individual friction sensitivity, surface exchange, reciprocal force and work,
original connections, numerical-accuracy refusal and exact retry.

Its `apply_source_refresh` adopts the original circular-line material endpoint
only against the complete current predecessor. Both physical footprints, body
coordinates and velocities, connections, gap column, time and historical
observations remain intact. The original source error estimate and an independent
combined mechanical/surface parameter-work balance must fit the native energy
gate before source, state and load line publish together. Subsequent native
motion uses the new surface tension, wetting thresholds and line resistances.
The returned signed surface parameter work is a separate account; it is not
automatically deposited heat or physical contact-line loss.

## Constitutive graphs on the finite thermal timeline

`constitutive_graph::enthalpy::step_candidate` closes public material-law feedback on the existing `DepositionSimulation` and `BackwardEulerStep`. Each endpoint trial binds actual caloric temperatures and the original physical interval duration to declared graph inputs, solves the original dimensioned material equations using their consistent tangents, and deposits explicitly bound watt-valued graph outputs once into the receiving cells. An outer energy-coordinate Newton solve includes the converged material response in its numerical derivative. All trials borrow the same committed graph history; the returned joint candidate carries the thermal endpoint and proposed law history together. Fresh local and total joule gates and an independent external-supply evaluation precede return. Original material-source failures, events, cancellation, and aggregate numerical work limits preserve both committed histories.

`PreparedHeatNetwork::graph_step_candidate` composes this feedback with the existing material heat-network resolution, including source-integrated conductivity attachments, interfaces, and radiation. Its thermal-only admission refuses electrical circuit bindings whose losses would otherwise be omitted or counted twice. Closed-enclosure radiative roundoff is separately reported and excluded from the independent external-supply account. The adapter retains fixed mass, geometry, and topology within the interval, original caloric and latent storage, and first-order backward-Euler rates; it does not certify passivity of arbitrary authored graph laws, supply phase kinetics, or remove the need for independent time-step refinement.

## Continuing ambient-fluid changes on native rod trajectories

`RodMotion::rebind_submerged_fluid_candidate` updates an already bound material-fluid environment at the accepted clock without changing pose, linear momentum, centroid spin, mass, inertia, supports, or step count. The candidate reevaluates the native response at the actual moving pose. `RodFluidRefresh` reports both coherent transport snapshots and the change in hydrostatic potential relative to the original reference geometry; that change is external parameter work, with zero instantaneous impulse and no invented heat. It does not represent the thermodynamic cost of replacing a finite reservoir. Subsequent native kicks include gravity and buoyancy forces and their actual rotating attachment couples in the external linear and angular impulse accounts; conservative potential work is not counted again as applied-wrench work. Failed admission or cancellation returns no changed trajectory. Focused environment tests cover data-driven sinking, neutral buoyancy and floating; analytic piecewise acceleration across a density change; eccentric strands with different mass and displacement centroids; gas-pressure changes; stationary reactions; and rejected or cancelled refreshes.


## Simultaneous adhesive events on complete elastic networks

`contact::adhesive_events::linear_network::multiple` resolves up to eight
independent scalar junctions on the same complete mass-normalized network.
Every original mode, modal damping term, rank-one spring/damper connection and
rest extension remains in the shared matrix-free cell propagator. The scalar
owner uses that same physical propagator and resistance quadrature.

All sticking reactions solve one normalized Gram system; its orthogonal
projection constrains only the active contact coordinates. Sliding uses each
held kinetic capacity. Dense-polynomial roots locate the earliest stop or
release across every junction, including hidden intervening stops. Events
within the declared time resolution are applied together. Zero-velocity
continuation retains the largest feasible sticking set; equally sized
admissible alternatives refuse, rather than depending on contact ordering.
Dependent or unresolved contact columns are refused without regularization.

Per-junction slip, friction loss, impulse and impulse time moment remain
separate from original modal/link damping and external work. Work uncertainty
includes event timing, constraint projection, propagation and normalized-Gram
conditioning; it must fit the original energy gate. These are numerical
estimates, not interval event certificates. Active-set factors are cached,
workspaces are preallocated, and aggregate contact, cell, operator, root and
phase-trial budgets bound computation. Accepted state and phases belong to the
outer physical transaction; refusal and cancellation preserve identical retry.

The lane holds source capacities, local work columns, coefficients and external
loads during each step. It does not supply vector Coulomb cones, redundant
contact force distributions, changing normal geometry, nonlinear material
memory, or acoustic back-reaction. G1/G3/G4 tests exercise analytic nonorthogonal
stops and release, independent dense-ODE work, simultaneous/hidden events,
contact-unit and time-partition covariance, conditioning and ambiguity refusal,
limits, cancellation and deterministic retry.


## Simultaneous sourced adhesive events in the original network

`contact::adhesive_network_events::multiple::MultipleAdhesiveNetworkEventModalSystem`
retains up to eight independently sourced scalar junctions on one complete
original modal network. Each site supplies its own ordered normal/shear source,
held indentation and complete work column. Original DMT real area and junction
strength determine its static and kinetic capacities. Coupled event propagation
selects holding reactions, stops and releases jointly; contacts are not advanced
as separate pairwise mechanical steps. Per-contact phase, accumulated slip,
friction loss, impulse and time moment remain available alongside global work.

Single and multiple owners share one original native-state transaction. It
restores every original component under its budgets, independently evaluates
pressure and storage, retains every spring/damper/rest extension, and computes
each original connection's force and work port. All individual source, slip and
interval-speed limits and the original energy budget, including propagated
numerical work uncertainty, must pass before any accepted body, phase, frame or
clock changes. Late original-body or energy refusal and cancellation preserve
all junction histories for identical retry.

The fixed scalar work columns must be independent under the admitted normalized
Gram pivot. Dependent constraints are refused without regularization or silent
merging. Two scalar contacts are not a vector Coulomb disk. Normal geometry,
source conditions and external force are held during a step; evolving frames,
normal contact, nonlinear memories and radiation feedback remain separate
physical owners. G1/G3/G4 tests include simultaneous three-body stopping and
restart, constrained internal elastic motion, per-site order covariance,
zero-strength recovery of the original complete prepared network, reduction
to the single-junction owner, source/rank admission and joint rollback.


## Finite-rotation native rigid cohesive patches

`cohesive_contact::rigid::RigidCohesiveSystem` evolves two original
`fs_mbd::RigidBodyState` owners and their original full `MassProperties` through
a single finite cohesive patch. Every site retains its source-bound original
point, irreversible history, reference area, body attachment arms and material
frame fixed in body A. Ordered receiving material identities and actual initial
geometry are mandatory; construction does not infer or reset a bond.

A joint twelve-variable wrench solve replays both original native Cayley steps
and all original objective constitutive trials. Its fixed dense dimension does
not grow with quadrature; material work is linear in the site count per residual
evaluation, with preallocated site scratch, prepared LU and explicit evaluation
budgets. The actual midpoint arms transmit patch moments. The original frame
couple acts on body A, and spatial moments use the inverse arithmetic rotation
average, tau_body = R0^T (M - c cross M / 2). These conventions retain reciprocal
force, total angular impulse and the original finite-rotation work identity.

Actual native storage, original cohesive storage, intrinsic fracture loss, bond
viscosity, external force/couple work and gravity are checked independently.
Body states, every point, accepted observations and the clock publish jointly
after force/moment, energy, momentum, source, geometry and final cancellation
gates. Late refusals leave an identical retry possible. Failed bonds retain
fragment motion and damage; the original compression penalty provides
frictionless reclosure within its declared domain.

Ordinary steps keep source conditions fixed; explicit temperature updates are
described below. Attachments and topology remain fixed. This is a bounded local
branch solve over straight local constitutive paths, not a proof of global
softening-branch uniqueness or intrastep trajectory accuracy. Rotation-chart
and original material/rate limits remain mandatory. Bond creation, healing and
collision discovery remain outside this owner. Optional damaged-area friction
uses the joint original-body owner described below.
`tests/rigid_cohesive.rs` checks an independent axial oscillator, global frame
covariance, finite-patch peeling with separate fracture/viscous losses, full
native free motion after failure, reclosure, gravity/external work, source and
geometry refusal, work-budget/cancellation refusal, strict energy refusal, and
late body-budget rollback with bit-identical retry.

### Temperature-dependent adhesion on the original native rigid pair

`RigidCohesiveSystem::step_with_temperature[_under_gate]` accepts one original
source-temperature successor per existing site, in site order. The original
objective owner transfers common- or mixed-mode history once at fixed accepted
geometry before any Newton evaluation. Reusable per-site scratch supplies these
same transferred points to the existing twelve-unknown native wrench solve;
mechanical probes do not repeat thermal damage or coefficient work. Original
source claims, ordered materials, physical mode, site areas, attachment arms,
full body inertias, poses, momenta and accepted clock remain bound.

The final original body audit balances kinetic/gravitational plus cohesive
storage change and all physical losses against mechanical load work and signed
`parameter_work_j`. `fracture_dissipation_j` includes thermal and mechanical
fracture exactly once; `thermal_damage_dissipation_j` identifies its thermal
subset. `thermal_stored_change_j` is the fixed-jump change before either body
moves. Per-site observations retain the same distinctions. Parameter work is
not a mechanical port or an inferred heat source. Each initial history transfer
consumes one original material-evaluation allowance. Following ordinary steps
retain accepted source temperatures with zero new thermal work.

Both bodies, every site/source, observations and clock publish only after the
existing force, rotation, geometry, energy, momentum and final cancellation
checks. Partial source staging and late failures leave accepted state unchanged;
the thermal scratch selector resets before any Result returns. Enabled
damaged-area friction requires the explicit paired-source operation below;
the cohesive-only temperature API never silently freezes its coefficients.

This is a prescribed-temperature first-order staggered composition. Thermal
expansion, heat transport, fracture-heat partition, cure, measured thermal
calibration, continuous temperature-path validation and global softening-root
uniqueness remain outside it. Focused native regressions compare independent
heterogeneous-temperature midpoint oscillators, finite rotating mixed-mode
patch work/momentum and frame covariance, source/late-body/cancellation rollback,
material budgets and the friction boundary. Their runtime execution remains
blocked here by the unchanged pinned executor dependency's memory requirement.

`step_with_temperature_and_friction[_under_gate]` accepts one cohesive and one
residual-friction successor for each original site. `ResidualCoulombSource`
refresh pins its original card, ordered surfaces, complete selected claims,
absolute-temperature convention and held coordinates, and admits the updated
static/kinetic pair together with `0 <= mu_k <= mu_s`. Every paired cohesive and
friction query must agree exactly. No coefficient is inferred from damage,
cohesive strength or temperature alone.

The existing rigid-Coulomb specialization stores no elastic friction energy.
Consequently changing mu at fixed jump supplies zero friction parameter work
and zero friction loss. The original cohesive transfer alone accounts for
changing compression/bond stiffness and any additional fixed-jump damage;
compression storage and fracture release are never duplicated in friction.
The same transferred material point, retained slip coordinates, prior fracture
history and static/kinetic occupancy enter every subsequent joint Newton trial.
Actual moving compressed endpoints retain kinetic capture, including newly
damaged material. The accepted law, source reference and state continue after
publication; no pristine friction or bond is constructed during thermal refresh.

Preallocated friction-temperature scratch shares the existing transaction and
material-work budget. Exactly one constitutive temperature transfer per site
precedes the original simultaneous body/friction solve; each final local audit
uses that same point and the new coefficients. Source mismatch, invalid static/
kinetic data, exhausted final-audit budget and cancellation leave original
bodies, all sources/histories, observations and clock unchanged. Mechanical
friction work remains the sole friction-loss account; the prescribed-temperature
extension adds neither a thermal recipient nor an entropy/heat-partition model.
Focused source/native regressions cover static heating/cooling with zero heat,
analytical two-tangent sliding with separate coefficient work and friction loss,
nonzero retained damage/slip, source identity and rollback with ordinary/thermal
retries. Native runtime execution retains the pinned executor memory limitation.

## Molecular gas refresh on continuing DG ports

`DgGasRefresh::UniformMolecularDiffusion` selects the native uniform-background molecular/thermal entropy-and-compression transfer in the existing `prepare_dg_gas` transaction. It supports the original exact-mass and spatial wall routes, including the explicit spatial quadrature used by `DgSourceInput` gas fields; it never flattens spatial histories. The original molecular source law, species fractions, geometry and reservoir support remain bound. The same option is forwarded by existing one- and two-sided Bernoulli/DG gas-refresh APIs.

The complete port transaction rebuilds the native propagator and algebraic port relations while preserving independently sourced branches, terminal channels and every original exterior MIMO pole state. Existing equilibrium-pressure and total absolute parameter-work checks still apply. Interior material changes do not reidentify the exterior medium. Parameter work is charged once, separately from molecular/thermal dissipation; this operation does not itself solve background heating or species transport. A regression in `dg_port::gas::tests::molecular` exercises nonzero molecular histories with an independently sourced branch/channel and full mutual rational terminal. Its full-crate integration execution remains separate from the focused native DG constitutive/continuation tests.

### Source-material walls in molecular pressure feedback

`MolecularDrivingPoint::new_with_material_wall` retains one
`Arc<LocalMaterialWall>` and passes that same immutable source law into
`MolecularDuct::new_with_compliant_wall`. The material identity, thickness,
axial constraint, supported anisotropy and elastic/viscous/Maxwell response
remain with the existing wall owner. Feedback, terminal flow and all finite
chimney observations therefore consume the same jointly loaded gas operator.
`material_wall()` exposes the original binding for inspection.

The optional gas shear/thermal resolution remains independent of mechanical
wall displacement. The existing characteristic FIR ladder and pressure
reflection, finite-hole loads, opening fractions, wavefront factoring and
Rayleigh outlet observation are unchanged. The positive material surface
mass makes its high-frequency admittance decay as 1/omega, so the leading
wavefront remains the original frozen molecular gas wavefront.

Every actual radius and queried frequency must satisfy the material wall's
domain. The existing near-DC probe also requires material support; missing
data cannot fall back to a rigid wall. The caller still declares the exterior
gas independently. This is a single local membrane material over the admitted
duct and chimneys, without nonlocal shell coupling, exterior wall radiation,
nonlinear gas evolution or a new performance-file grammar.

The focused integration tests independently assemble a Debye caloric law and
a membrane-loaded cylinder, compare pressure reflection and physical outlet
flow, vary gas and solid relaxation independently in the receiver response,
and check common source-band refusals. These tests are retained for the
canonical coupling build; the session's executed native slice verifies the
underlying molecular/compliant-wall operator but does not execute these
coupling/FIR tests.

## Finite source-material cavity walls

`SharedCavity::new_with_cylindrical_material_wall` couples an explicit cylindrical gas region to a finite solid wall using the original `PreparedThermalTransport` source, evaluated at the reference gas temperature. Radial finite-volume capacities and Fourier conductances include gas, wall, and any original local molecular reservoirs in one caloric operator before pressure elimination. Its positive compression realization joins the existing exact mechanical runtime and harmonic response, with one accepted state, pressure, availability, dissipation, and rollback account. Gas and wall cell counts are independent spatial refinement axes; the cold realization is limited to 128 total thermal coordinates and rejects unresolved or unrepresentable spectra. The exterior wall is held at the reference temperature and endcaps are insulated. The scope is small perturbations, uniform acoustic pressure, fixed geometry, and frozen reference coefficients; this path does not infer mean heating, phase evolution, thermal expansion, interface resistance, exterior convection, or spatial accuracy from the declared frequency band.

### Complete vibrating bodies on the reed-driven ODE duct

The reed ODE renderer consumes the original `PlateBankImage`, retaining compact
and Vk modes together, full material damping and every accepted internal state.
The signed surface-volume port supplies reciprocal pressure reaction independently
of the external point-force port, which remains unforced. The complete original
plate observer uses the same state prefix. The existing inlet-pressure plus plate
observation convention remains; this change does not supply a new exterior outlet
microphone model.

For an exactly quadratic body/duct, a quasistatic reed uses the same complete
`PreparedExactZoh` held-volume-flow owner as the massive reed. The original strict
aperture root closes against the SAME interval's affine mean-pressure response,
with no delayed-pressure surrogate, hidden 5 Pa primer, extra initialization
step or render-length-dependent model selection. The original static opening
clamp and flow dead zone remain. Its physical flow is reconstructed from that
constitutive law; positive finite mean body impedance and the actual accepted
body pressure/work must agree before publication. Complete state/output buffers
are reused.

This connection gives exactly `body_work = dt * mean_pressure * held_flow`.
The massless valve's independently evaluated jet loss and mouth/body work obey
`jet_heat + body_work = mouth_work`; zero mechanical storage remains the declared
quasistatic model. The prior held-wave Robin connection included a within-window
pressure variance term absent from the valve's mean-flow work account. The
held-flow connection removes that coupling loss without modifying original
material/acoustic resistance or smoothing the aperture law. It does not certify
the physical adequacy of a massless model near a reed resonance.

For a massive reed and exactly quadratic body, the prepared complete body now
supplies its current-window held-volume-flow response `p_bar = free + B*q`.
`reed_bore::loaded::PreparedReedLoad` closes the original moving aperture, swept
flow and lay reaction against that response through the existing scalar contact
solver. The original linear reed mechanics evolve exactly under the interval's
held pressure/contact force; the aperture uses the midpoint gap. The body's
actual accepted mean pressure and work must agree with the condensed junction
before both candidate states are published. No previous-sample pressure drives
this quadratic path. Nonlinear bodies use the complete exponential
discrete-gradient path specified below; they are not inferred quadratic from
probes. Existing dense-state limits remain. Source histories seed a private
one-shot render snapshot.

The loaded reed retains original contact potential, exact integrated viscous
heat, nonadhesive contact loss, and Bernoulli jet loss as separate quantities.
It checks `deltaH_reed + viscous + contact + jet + body_work = mouth_work`.
The exact held-flow body's original storage, material/acoustic loss and port work
close the other side of that same interval. A finite strictly positive mean
body impedance is required; there is no impedance clipping, state projection
or substituted numerical damping. The source scenario's existing effective
area and lay calibration are retained, without implying measured specimen
authority. Exact linear propagation does not imply exact nonlinear continuous
motion or root uniqueness; the nonlinear aperture/contact step needs temporal
refinement.

Preparation owns the existing exact maps and contact storage. Warm reed trials
reuse the existing scalar solve and allocation-free original contact response;
exact blade velocity and viscous-heat contractions no longer allocate temporary
scaled states. A 1,000-candidate probe spanning open and contacting states
observes zero allocations. This is an allocation result, not a real-time timing
qualification.

Focused native tests execute the complete renderer with original gas,
viscothermal Foster, modal material, contact and numerical owners. They check both
reed selections with mixed modes and non-diagonally-dominant damping, independent
drive/surface maps, zero-blow silence and warm-memory response, render-length prefix
identity, same-window aperture pressure/flow against independent physical RK4,
and preservation of the original direct Gonzalez clock refusal. These are
numerical/structural checks, not measured-instrument, complete workspace or
real-time qualification.

## Simultaneous Hertz surface and material rolling bulk heat islands — v0

`friction_heat::hertz_rolling::multiple::propose_rigid_hertz_rolling_island_heat`
joins one immutable prepared native five-coordinate contact island to the
original finite enthalpy network. Each original surface-friction source is
resolved at its own solved common-interface temperature and explicit pressure.
Each original constant-coefficient rolling material card is admitted at its
deforming sphere's actual bulk temperature. Interface identity must agree
between the original rolling and surface-friction owners.

Every patch routes its surface dissipation once through its distinct original
thermal link. Bulk losses from multiple patches on the same sphere are summed
into one original sphere cell, retaining one heat capacity. Different spheres
retain distinct bulk temperatures, inventories and heat destinations. Admission
requires the native whole-sphere mass, declared body/material identities, the
correct original route side and a caloric centre coincident with the native COM.
Conflicting mappings, a shared inventory for different spheres, duplicate
active sphere cells and duplicate thermal links refuse. Canonical permutations
between original patch IDs, thermal links and unique sphere bodies are
prepared once and retained throughout the nonlinear solve.

The original thermal network retains ordinary conduction and authored outside
heat. Signed outside heat is measured before injecting internal bulk rolling
loss, preserving its physical boundary and small signed contributions. Every
trial starts from unchanged accepted native momenta, original rolling histories
and caloric state over the same actual representable interval. The callback
supplies original endpoint conductivities and other heat only; interface work
and bulk work have distinct original owners and destinations.

Acceptance independently checks actual
`Delta H - Q_outside - sum(D_surface) - sum(D_bulk)`,
`Delta K + Delta H - Q_outside`, and the original native/network residual sum.
Only a complete mechanical/history/source/thermal candidate is returned.
Source failures, late audits, exhausted trials and cancellation leave accepted
owners unchanged and permit exact retry. Original cell/unknown, trial,
nonlinear, mechanical and physical-residual budgets remain explicit.

This retains frozen normal loads and poses, independent nonoverlapping local
sphere-deformation fields, the original rolling-card applicability box, one
uniform caloric sphere cell, and common-face backward-Euler transport.
Temperature admission does not create temperature-dependent viscosity.
It adds no spatial thermal validity, retained partial-slip memory, entropy
classification of injected bulk work, root-uniqueness or experimental claim.

Six focused G1/G3/G4 groups in
`tests/articulated_friction_heat/hertz_rolling_multiple.rs` use independent
coupled torque equations and a two-face thermal root, check one shared sphere
inventory, Hertz-induced rolling from rest, distinct original source
temperatures, signed outside heat and exact patch permutations. They also
check stationary conduction, zero viscosity, invalid inventories/interfaces,
late source failure/cancellation/budget retry, and two distinct deforming
spheres with nontrivial native/patch/link ordering and independent per-cell
heat destinations. All 23 existing groups in the original articulated-friction
heat test target continue to pass, for 29 groups in that target.

### Source gas in spatial thermoviscous interiors

`cavity_modal::interior::thermoviscous_interior` binds the original spatial
Helmholtz pressure solver and its viscous/thermal wall form to one coherent
`GasState`. Density and sound speed must match the FE volume medium bitwise;
viscosity, thermal conductivity, isobaric heat capacity and gamma are copied
from that same material evaluation. No dry-air constants or alternate gas law
are substituted. The original FE owner retains geometry, reduction, frequency
band and boundary-work semantics, including raw wall power up to roundoff.
This is a quiescent homogeneous frequency-domain interior with stationary
isothermal thin boundary layers. It does not add molecular relaxation,
finite-conductivity walls or narrow-gap linearized Navier–Stokes.

Six focused new native cases passed against actual FE/sparse/modal/material
source owners: four boundary-layer cases and two source-gas integration cases.
The complete Cargo workspace was not built in this environment.

### Spatial molecular gas changes in continuing port networks

`DgGasRefresh::SpatialMolecularDiffusion` routes the existing main/selected-branch
transaction through `DgDuct::rebind_spatial_molecular_gas_field`. The original
DG owner preserves local thermal/molecular entropy and weak compression/momentum
on unchanged spatial quadrature, while the complete port transaction retains
branch/rational-terminal histories and charges the aggregate absolute parameter
work. Existing common-pressure junction/exterior checks and all-or-nothing
publication remain in force. The background may change temperature and pressure
but must retain station positions and complete species fractions. This adds no
transport law, interpolation of histories or new exterior gas assumption.

Four new native physical-transfer tests and five existing uniform-transfer
regressions passed using actual DG/material owners. The added complete branch,
channel and mutual-terminal integration test is present but was not executed
because the full coupling dependency graph is unavailable in this checkout.

### Continuous molecular fields in feedback and opening sound

`MolecularDrivingPoint::from_gas_field` and its material-wall counterpart bind
the original primitive `AxialGasField` and molecular law to the same characteristic
pressure feedback, terminal flow and finite-chimney mouth observations. Axial
Magnus resolution, radial diffusion resolution and FIR length remain separate
axes. Inlet loading uses the exact inlet gas, exterior radiation uses the explicit
exterior gas, and terminal/side-mouth causal delays use the core's frozen-slowness
quadrature and each chimney's local gas. Uniform field arithmetic is preserved.
The source-resolved local material wall is shared with the propagation owner;
its positive surface mass preserves the gas-frozen wavefront.

The existing whole-field coverage, source band, common pressure, internal mutual
hole and work-cap refusals remain. The existing reed constructor accepts this
prepared driving point. The field is prescribed and frozen during playback;
this does not solve species transport or add nonlocal structural/exterior loading.
The two core wavefront regressions ran natively. Three added full feedback/FIR
consumer regressions are present but unexecuted because the complete coupling
and material-wall dependency graph is unavailable in this environment.

### Independent massive reed/body coupling checks

The actual full Foster-bore and massive-reed equations are compared with an
independent fine RK4 reference. Relative interval-mean pressure errors at
24/48/96 kHz are 2.215118e-4, 5.541245e-5 and 1.385526e-5, approximately fourfold
reduction per refinement for the smooth constant-pressure fixture. Fourfold
reed mass at unchanged static compliance changes actual bore pressure by 2.14463%
under that test's relative norm. Independent analytical damped-oscillator and
positive Simpson heat references cover contact entry, penetration and release,
with original contact potential, jet heat, mouth/body work and refusal replay.
Existing contact regressions remain active. These bounded numerical checks do
not establish measured-instrument fidelity or arbitrary-clock nonlinear accuracy.

### Molecular gas fields in authored reed performances

`ReedPerformance` accepts the existing `molecular`, `molecular_resolution` and `molecular_wall` records after an explicitly continuous `gas_field` or `gas_steady` block, or with a field injected through `from_bytes_with_field_resolution`. The parser retains the complete `RelaxationAbsorption` declaration and admits it against the actual field inlet. The continuous molecular owner then reevaluates the same law at every original field/geometry quadrature point, including species absent at the inlet but present downstream. The legacy homogeneous path retains its original evaluated-state constructor and arithmetic.

The exact field inlet supplies excitation and feedback; the exact outlet supplies the default baffled exterior, with explicit exterior records retaining precedence. The same molecular operator drives the ordinary reed voice, terminal/chimney volume flows and pressure observation. Authored-segment inlet samples remain read-only metadata; they do not replace continuous quadrature. Field/source/steady-solve receipts and observer-law selection retain their existing meanings. Injected fields retain their complete source owner while `gas_records` continues to mean authored records.

Molecular fields require `GasFieldResolution::Continuous`; midpoint input refuses instead of being silently reinterpreted. Explicit segment gases remain unsupported by this molecular grammar. Existing source domains, composition limits, common-pressure admission, finite work caps, internal mutual-hole refusals, FIR approximation and frozen-field scope remain in force. No new grammar, transport model, source fitting or audio effect is introduced. The new `molecular_input::field_tests` cover actual authored/injected field-to-pressure/block replay, law sensitivity, downstream source/domain refusals and the existing stationary thermal/species generator feeding the molecular performance. These full coupling/render tests require the complete workspace test lane.

All three new full coupling/render regressions are unexecuted in this environment.
Rust formatting/parser checks passed; the constituent molecular field and wavefront
operators have separate native evidence, without claiming a full render run.

### Original work through a quasistatic valve and a warm complete body

The first-window reference integrates the original body J/R/G under constant
physical volume flow independently of the prepared maps. A mixed compact/Vk body
with off-diagonal damping and warm Maxwell history is checked at zero and 100 Pa
mouth pressure by independently integrating its full states, pressure and original
R heat. Both body storage plus heat equals bore work and body storage plus heat
plus jet loss equals mouth work. At zero mouth pressure the original body motion
drives reverse valve flow and pays for both loss channels. Direct helper checks
also cover signed reverse flow and complete static closure. Existing massive
contact/work, refinement, state-prefix and nonlinear-refusal checks still pass.


## Joint two-tangent adhesive friction on complete linear networks

`adhesive_events::linear_network::vector::VectorLinearNetworkFlow` advances
the original mass-normalized modal coordinates and every original passive
spring/dashpot link with simultaneous Euclidean Coulomb disks. Each contact
has two independent tangent columns and one capacity; the monotone rung has
equal static and kinetic strengths. It has genuine set-valued sticking and
direction-changing vector slip, without per-axis force clips or smoothing.

The complete backward-Euler operator is condensed to at most sixteen contact
unknowns. A prepared mechanical factor and mobility are reused at a fixed cell
size; analytic disk projection derivatives drive a bounded semismooth Newton
solve, with a projected-gradient fallback. Workspace allocation occurs at
admission. Unit-column Gram pivots require independent, conditioned contact
coordinates; no dependent reaction choice is assigned physical authority.

Each cell measures original modal/link damping, contact impulse and its first
time moment, per-contact slip and friction work, and original stored energy.
The exact discrete loss `0.5*|delta_v|^2 + 0.5*delta_q^T*K*delta_q` is reported
separately as numerical dissipation, including each connection's spring share.
It is never material damping, friction heat or physical fracture loss. An
independent uncondensed equation check and Coulomb work gap contribute to the
numerical work estimate. Explicit cell, nonlinear, cancellation, physical-state,
energy-accuracy and cumulative numerical-loss gates precede return. Failed calls
retain no accepted history and retry deterministically from caller-owned state.

This is first-order discrete mechanics; cell refinement remains necessary.
Work estimates are not interval or trajectory-error certificates, and link-force
maxima cover discrete endpoints. No exact event time, unequal static/kinetic
strength, evolving normal geometry, moving tangent frame, material memory or
nonlinear network is supplied. Nine G1/G3/G4 tests cover analytic oblique stops,
disk-versus-box behavior, joint momentum/sticking, independent full-network RK4
and first-order convergence of state and numerical loss, basis/contact-order
covariance, internal elastic motion, partitioned impulses and moments, and
strict physical/numerical/cancellation refusal with bit-identical retry.


## Original native mechanics with sourced vector adhesive friction

`adhesive_network_events::vector::VectorAdhesiveNetworkModalSystem` binds every
two-tangent disk to its original ordered DMT/shear sources. Real contact area
and the original junction strength determine capacity. This owner explicitly
requires equal static and kinetic source strengths; it does not flatten a
velocity-weakening law into a monotone one. All native modes, modal damping,
Kelvin-Voigt links, rest extensions, pressure observations and physical limits
remain active. Each contact must supply both complete independent work columns.

The vector kernel and scalar event owners share original native state staging.
Actual component restoration, independent storage/work checks, each original
connection's force/work, source/speed/slip limits and propagated numerical work
uncertainty precede joint publication. The vector owner's backward-Euler loss
has its own `numerical_dissipation_j` field in both native and adhesive frames,
and a separately named numerical-loss entry in the native window audit. Physical
modal damping, link damping, friction, external work and numerical loss remain
distinct; `CoupledModalFrame::total_dissipation_j` contains physical losses only.
Each connection's discrete spring loss also enters its actual work transfer.
No numerical loss is assigned to heat, and no residual is erased as a loss.

Eight G1/G3/G4 tests cover diagonal loading outside the disk, joint sticking
with continuing internal elastic motion, tangent-chart rotation and contact
permutation, individual source causality, first-order refinement toward the
original complete prepared network, separate numerical/work-accuracy limits,
source and rank admission, cancellation and late original-body refusal with
bit-identical retry. The 26 existing scalar native event tests also pass.
This lane retains fixed normal indentation, fixed tangent geometry and source
conditions, with bounded internal timesteps. It does not add rotating frames,
normal impacts, velocity weakening, nonlinear constitutive memory, dependent
reaction selection or exact continuous event times.


## Simultaneous circular wetting bridges on original native mechanics

`capillary::multiple::MultipleCircularWettingModalSystem` binds up to 32 original
asymmetric circular contact-line bridges to one original modal network. Every
bridge retains its two independently sourced surfaces, physical volume,
footprints, pinning/sliding state and loss channels. Complete reciprocal gap
columns condense the full native displacement response, including all
cross-bridge terms, original component modes and bilateral connections.

A bounded simultaneous Newton solve uses the original analytic force tangents,
prepared dense LU and backtracking. Probes always start from accepted footprint
states over the same physical timestep. One aggregate source-evaluation ceiling
includes rejected probes and reserves the final actual-motion replay of every
bridge. All bridge/Newton/candidate buffers are allocated at admission; the
original native propagator retains its own allocation behavior. Inverse-Jacobian
gain screens ill-conditioned local continuation and propagates original force
uncertainty through every cross-bridge compliance. These remain numerical
estimates, not interval, global uniqueness or stability certificates.

One native stage applies the sum of solved forces. Every original bridge is
reevaluated at its actual staged gap, and body motion, all footprint histories,
gap coordinates, diagnostics and sample clock publish together after source,
force, compatibility, complete work/energy and cancellation checks. Original
surface storage, native loss, per-line viscosity/pinning and external/prescribed
work remain distinct. No heat recipient is invented.

Seven focused G1/G3/G4 tests exercise analytical three-body impulse/motion,
coupled two-bridge work, independently pinned footprints, six original elastic
modes with pressure and connection loss, reduction to the existing scalar owner,
bridge permutation, cross-bridge source causality, aggregate work/conditioning
limits, source identity, cancellation and strict/late-native atomic retry. The
18 existing circular and planar/axisymmetric native tests also pass. Each bridge
retains the original connected, acute-angle, nonpositive-pressure, quasistatic
liquid scope. This owner does not model liquid exchange/coalescence between
bridges, changing contact geometry, bridge rupture, liquid inertia, evaporation,
thermal deposition or exact within-sample trajectories.

`MultipleCircularWettingModalSystem::apply_source_refreshes` accepts a bounded
nonempty subset of those same original circular-line endpoint candidates in one
atomic transfer. Every selected bridge must match its complete current history;
duplicate indices refuse. Canonical bridge-index reduction retains input-order
independence, all cross-bridge mechanical compliance and untouched bridge states.
A late predecessor, source, energy or cancellation refusal publishes nothing.
Native regressions cover thermal unpinning during continuing motion, combined
surface parameter work, batch permutation, local updates and unchanged body,
footprint, clock and observation history on both successful adoption and refusal.


## Damaged-area friction on original rotating rigid cohesive patches

`cohesive_contact::rigid::RigidCohesiveSystem::enable_damaged_area_friction`
attaches one original `ResidualCoulombSource` to each existing bonded site,
without replacing either native body, changing its full inertia, resetting
damage/history or advancing the physical clock. Original ordered surface cards,
query conditions and material assignments must agree. The explicit constitutive
selection is the shared damaged-area rigid Coulomb mixture: each tangent disk
has capacity `mu * area * integral(damage * compression_pressure)`, evaluated
analytically on the same original material path. Tensile cohesive traction,
surviving viscosity and endpoint damage are not substitutes for that measure.

Two local tangent-force unknowns per site join the original twelve body-wrench
unknowns. Every disk, body pose and material point is solved simultaneously.
Static capture uses the shared law's static coefficient, or kinetic capture
when already sliding; a yielded capture solution resolves the entire patch at
its kinetic coefficients. Numerical cone scaling introduces no elastic spring.
The explicit `velocity_absolute_tolerance_m_s` is independent of timestep.
Actual objective velocities from the original native bodies select kinetic
capture and detect resolved endpoint breakaway, even when a short interval's
slip is below the displacement tolerance. The same fixed speed gate applies
after a source-temperature transfer. A static path requires both endpoint
velocities and net slip to be unresolved; resolved moving paths retain their
original kinetic `mu * damaged_compression_force * slip` loss under refinement.
After an actual stop, endpoint occupancy restores the static capture disk for
the following interval without changing the audited preceding kinetic loss.
Prepared damped least-squares Newton corrections admit redundant sticking sites
without changing physical compliance or the final force equations. Individual
reactions may be a numerical representative of an indeterminate distribution;
neither their uniqueness nor global softening-branch stability is claimed.

Friction force uses the original arithmetic rotation-average pullback and its
balancing frame couple. Actual native linear/angular momentum and kinetic plus
gravity work therefore use the same finite-rotation geometry as each original
cohesive point. The shared constitutive audit independently checks the actual
local displacement, full two-dimensional Coulomb disk and dissipation; a second
spatial audit checks its objective work. Fracture, surviving bond viscosity and
friction losses remain separate, and friction is not automatically called heat.
The original surface compression energy is counted once. Complete failure
retains physical fragments and residual friction on reclosure, with no healing
or invented bond creation.

Site count bounds dimension at 140. Dense Jacobian, factorization and site
scratch are allocated at admission; no per-step friction buffer allocation is
required. Physical evaluations are linear in site count, while the joint
numerical differentiation and dense solve have their explicit bounded costs.
All probes, branch switches and final material audits share the original work
ceilings and cancellation gate. Bodies, all bond/friction states, observations
and clock publish only after the complete original audit; a refusal preserves
all accepted state and deterministic retry.

Eight focused G1/G3/G4 tests cover static holding between the static and kinetic
thresholds, an analytical two-tangent sliding trajectory, partial-damage area
scaling, independent quadrature through progressive fracture, distinct fracture/
viscous/friction losses, redundant quadrature partition, rotating spread patches
with full inertia, source/cancellation/late-energy rollback, and open/zero-friction
reduction to the original native dynamics. A resolved incoming slip cannot be
silently reflected by a zero-displacement midpoint capture: that case refuses
with a step-subdivision request. The same refusal applies to resolved endpoint
reversals with nonzero net slip. Two further native regressions check kinetic
work below displacement resolution at original and transferred temperatures,
an independently predicted actual stop followed by static holding, and
nonzero-slip reversal with complete rollback. Exact stopping times, general intrastep reversals,
contact discovery, changing attachment topology, wear/heat deposition and
arbitrary moving fracture-surface material reselection remain outside this owner.

## Prescribed-input graphs and finite thermoelectric stores

A `MaterialGraphResidual` may have zero algebraic unknowns and zero equations when every graph input is prescribed. The original zero-dimensional Newton owner still performs source admission and the final fresh graph evaluation; both evaluations count against the work allowance, and the proposed graph history advances once from the same committed history. Empty residuals cannot suppress source, missing-input, budget or cancellation failures. The original graph-to-enthalpy adapter supports native Kelvin and Ohmic nodes on the same finite thermal cells: terminal Peltier and bulk Thomson heat are signed reversible transfers, Joule heat remains native nonnegative dissipation, and the original material heat network retains Fourier transport. Current-driven cooling and voltage-constrained generation share the same graph/enthalpy solver; there is no thermoelectric-specific time integrator or invented storage certificate for arbitrary graphs. Focused public integration tests compare both paths with the original circuit thermoelectric owner, an independent constant-property four-cell heat equation, and the native energy/entropy balances under material-only changes and continued thermal history.

## Continuing rod motion under prescribed body temperatures

`RodMotion::rebind_material_temperatures_candidate` applies an admitted local section-temperature field at the accepted clock, preserving pose, momenta, kinetic frames, supports and integration-step count. It reevaluates native forces at the actual moving pose and returns `RodThermalRefresh`, whose signed fixed-pose strain-storage change is work supplied through the prescribed thermal parameter. Failed source admission, actual-pose strain bounds or cancellation return no changed trajectory. This operation consumes a supplied body-temperature field; it does not equate ambient reservoir temperature with body temperature or supply a finite thermal-reservoir balance. Sensible/latent storage, reciprocal thermoelastic heat exchange, phase evolution, cross-section expansion, mass transport and separately bound link heating remain with their respective physical owners.

### Sourced free expansion in native electrothermal circuits and field motion

`circuit_enthalpy::material::PreparedMaterialNetwork::with_free_expansion`
attaches original `PreparedFreeExpansion` paths to named whole-cell resistor
branches. `ExpandingMaterialNetwork` reuses the original material queries and
source pins, the SAME `MaterialOhmicNode` geometry/source admission as the
constitutive graph, and the native Kirchhoff/finite-enthalpy residual. At the
actual trial temperature, `L=L_ref*s_parallel` and
`A=A_ref*s_perp1*s_perp2`; `ResolvedConductor` retains resistance, current and
Joule power ownership. Conductivity is evaluated before its reciprocal, with
the original complete bundle intact. The new geometry lowering performs no
additional material query, and the source strain lines are prepared once.

The adapter's `step_candidate` preserves native charges/fluxes and routes each
branch's actual Joule heat once to its original caloric cell. Its
`field_motion_candidate` uses the same reciprocal electric/magnetic gap owner:
changed resistor geometry changes current, field force, mechanical work and
finite heating in one coupled solve. Resistor expansion does not change the
actuator's geometry or invent a force coefficient. Initial circuit sources and
geometry must match the actual accepted thermal state. Failed source queries,
invalid boundaries, nonlinear refusal or cancellation return no partial state.

This lane admits quasistatic, stress-free, wholly solid expansion at fixed
delivered mass. Reference thermal geometry and caloric charts remain intact;
changing density never changes mass. Expanded cells refuse occupied interbody
thermal contacts because those evolving interface areas and layer distances
require their own heat-transfer closure. Robin/reservoir exchange is admitted
only through the explicit source surface partitions below. Independent applied
watts remain allowed. Constrained thermoelastic work, arbitrary curved or
shearing surfaces, phase changes, skin effects and changing magnetic storage geometry
are outside this admission. Coupled temporal order remains the original first
order; no performance benchmark or unconditional passivity claim is made.

Four focused G1/G3/G4 groups in `tests/circuit_material_network/expansion.rs`
cover both electrical bases, all principal current directions, voltage/current
heating against independently solved endpoint equations, finite-inductor
current and energy, actual magnetic gap motion under strain-data-only changes,
source/boundary refusal and late cancellation/replay. Shared material lowering
also has a source-bundle parity/refusal regression. Native runtime execution
must be reported separately from parsing and test-target compilation.


## Physical Hertz stop/restart with original surface and bulk heat — v0

`friction_heat::hertz_rolling::collinear::propose_rigid_hertz_rolling_collinear_heat`
composes the native finite-body collinear stop/restart trajectory with the
original source-bound interface and finite thermal network. One immutable
native geometry is reused for all thermal root trials. Each complete trial
restarts the accepted momenta and opaque stop/material history, queries the
original Hertz source at its own common face temperature and pressure, and
admits the original rolling card at the sphere's actual bulk temperature.

The native trajectory's integrated surface loss enters the original massless
face once. Its independent integrated bulk loss enters the matched original
whole-sphere caloric cell once. The original mass, material/body identities,
COM, unique caloric inventory, geometry and work-route checks remain active.
Outside heat is measured before internal bulk injection. A true surface stop
can leave nonzero rolling loss, and ordinary thermal conduction continues for
the rest of the requested representable interval. No mean-rate or endpoint
power surrogate replaces the original integrated work.

Acceptance checks direct caloric closure
`deltaH-Qoutside-Dsurface-Dbulk`, direct native/caloric closure
`deltaK+deltaH-Qoutside`, and the independent original-owner residual sum.
The complete native stop continuation, original rolling checkpoint, source
receipt and caloric endpoint return together. A late source refusal, native
phase budget, thermal trial budget or cancellation returns no partial state.

For each trial, the mechanical coefficients are held at that trial's endpoint
temperatures and phase integration is exact for the admitted held model.
Storage and conduction retain the original backward-Euler window. Exact
stopping for those fixed coefficients is not an exact continuous-temperature
stopping claim; thermal refinement remains independent. The caller cuts source
or topology events. Geometry, normal load/radius and frequency remain held;
this owner adds no normal elastic work, general shear/spin or island event
trajectory, moving thermal footprint, variable viscosity or bulk entropy claim.

Six G1/G3/G4 groups in
`tests/articulated_friction_heat/hertz_rolling_collinear.rs` independently
reconstruct stopping time, scalar/exponential impulses and work, and the
finite two-capacity thermal root. They cover high static capacity, continued
post-stop conduction, temperature-dependent stopping, bulk rate variance,
post-stop rolling heat, authenticated continuation and late failure/retry.
The fixtures retain the unmodified original quasistatic material gates.

## Source-positive damping in the shared exact modal owner

The original exact modal owner retains every component and signed connection-dashpot source square before assembling its unchanged resistance matrix. Total resistance carries its checked native `ResistanceWitness` through material-memory and acoustic-load composition. Summed connection heat and each original link's work use the corresponding source factors against the same full augmented trajectory, while modal, material and radiation accounts retain their original owners. The source-factor path changes neither damping values nor mechanical dynamics and does not allocate a total heat remainder among channels. Original componentwise reconstruction, positivity, complete state/work admission and energy-balance checks still apply.

## Frame/fluid resistance: retained source squares and the reconstruction gate

`MaterialFluidFrame` builds its resistance as a sum of physical squares (the original material-damping ports, the low-frequency tail, one memory-slip square per spectral branch) and offers fs-phs the same squares as a `ResistanceWitness` (`PortHamiltonian::resistance_witness`, rows through `ResistanceWitness::rows`). The witness is attached only when fs-phs's componentwise reconstruction gate (`|gram - R_ij| <= 128 n eps |R_ij|`, with the local product sum as the scale at an exactly zero entry) admits those rows. When the gate refuses, the system is the unwitnessed one, `MaterialFluidFrame::resistance_witness_refusal` returns the gate's own refusal, and the gate is unchanged; an exact flow over such a system factors `R` itself and fs-phs decides whether that builds.

Two refusals are pinned by tests that repeat the gate on the retained rows and print the entries it objects to (`--nocapture`; `frame/harmonic/tests.rs` and `modal/relaxation/temperature/tests.rs`), measured on ts2 at tip 4d9153a5b:
- A dry four-mode viscous frame (`n = 8`, 6 source rows, so no QR compression): 3 failing entries, all between momentum coordinates. The first in gate order, `(4,5)`, has `R_ij = 1.0e-17` against a largest source term of `7.8e-4` (ratio `1.3e-14`); the worst, `(5,7)`, `R_ij = 7.6e-17` against `7.9e-2` (ratio `9.6e-16`). These off-diagonal entries sit 14 to 16 orders below the terms that make them up, so the gate's relative scale cannot be met; the numbers do not say why they cancel (orthogonality of the retained modes in the damping metric is the expected reason and was not tested separately). An exact ZOH over this unwitnessed system builds. The same frame with two fluid regions is admitted and carries its witness.
- A one-mode wet layered frame (`n = 10`, 12 source rows, QR-compressed): one failing entry, the structural zero `R(2,3)` between two memory coordinates. The source rows have no term in it; the compressed factor's Gram there is `5.1e-18` against an allowance of `3.0e-19`. That is fill-in of the compression at an exact zero, not cancellation of nonzero terms.

Frames whose rows the gate admits (the wet frames of the exact-flow tests in `modal/component`) keep the witness, and those tests assert it.

## Joint molecular and spatial cavity thermal response

`ReciprocalStringBank::with_material_shell_cavity_thermal_molecular` binds the original `RelaxationState` to the same pressure FE, heat eigenwindow and material-shell/string runtime. Each retained heat mode uses the shared positive caloric spectrum before pressure elimination; the positive unresolved pressure-space complement keeps the original adiabatic molecular spectrum. Their sum preserves frozen bulk storage without double-counting diffusion or dropping molecular reservoirs. Source identity, complete memory/local caloric work, native setup and instantaneous Nyquist gates remain binding. Pressure observation and the public normalized constitutive response use the same branch rows as mechanical force, storage and loss. Empty reservoirs recover the existing thermal trajectory exactly. This is a cold, frozen-mean, linear acoustic binding with isothermal solid boundaries and insulated mouths; it introduces no mean heating, phase/species evolution or finite solid-wall thermal store.

### Source-owned conductive walls in spatial acoustic interiors

`cavity_modal::interior::conductive_wall_transport` reads normal conductivity
and volumetric heat capacity from an original `ThermalTransportSample`, whose
prepared source identity remains on the supplied sample. It requires bitwise
agreement of solid and source-gas equilibrium temperatures and validates the
original gas coherence. No material-name constants, Cp/Cv conversion or second
heat-capacity model is introduced. The caller retains source preparation,
caloric constraint, direction and uncertainty ownership.

`thermally_conductive_interior` binds that sample, physical solid depth,
backside thermal condition and minimum tangential variation length to the
original `thermoviscous_interior`. The same FE medium, pressure basis, surface
matrices, response and input-work balance remain owners. The fs-feec normal
Fourier profile supplies gas/solid thermal boundary exchange at fixed mean
temperature; its geometry/admission and model-form limits remain in force.
This is not finite mean-body heating or a second application of gas-wall heat.
Two focused adapter regressions use original source cards/prepared transport
and source gases, checking material-only conductivity/capacity substitutions,
full FE work closure, and mismatched-temperature/incoherent-gas refusal.


## Nonlinear complete-body reed feedback

`acoustic_realize::reed_ode` uses `PreparedNonlinearZoh` for a body with original
nonquadratic storage. The complete quadratic mechanical, material-memory,
cavity, radiation and Foster-duct block evolves through the existing exact
linear owner. The original nonlinear remainder supplies a discrete gradient
closed implicitly against the accepted body endpoint; no mode or memory state
is dropped to admit the audio clock.

Each endpoint trial supplies its own same-window affine mean pressure to the
original moving-reed/contact or quasistatic aperture owner. Reed history stays
provisional throughout the body Newton solve. At the accepted endpoint the
original reed work must match actual body work and mean pressure, after which
the private source/body snapshot advances together. A failed observation or
numerical/source refusal publishes no source history. The original effective
reed face and reduced lay calibration are unchanged.

Body loss independently integrates the original resistance on the method's
held discrete effort. The energy account is for this exponential
discrete-gradient scheme, not an assertion of exact continuous nonlinear heat
or motion. Preparation and warm buffers are reused. The shared complete-state
capacity, bounded root work and positive mean impedance requirements remain;
contact and nonlinear temporal accuracy require refinement separately.

Both massless and massive public renderers pass at 48/96/192 kHz for a warm
mixed body with quartic cross-mode stress, full damping, Maxwell memory and a
complete Foster duct. Independent original-equation RK4 references show more
than 2.5-fold reduction in pressure error per clock doubling, with the finest
symmetric normalized L2 error below 0.002. Removing the original quartic force
changes pressure at the same initial state. A late-observer-refusal regression
preserves every warm source coordinate. The direct Gonzalez owner retains its
original clock refusal; this extension does not weaken that gate or claim
root uniqueness, measured-instrument accuracy or real-time throughput.

## Expanding-surface feedback in native electrothermal networks

`with_free_expansion_and_boundaries` binds disjoint principal-normal ambient
and reservoir face areas to that SAME expansion source, reference query and
frame. The original thermal surface and reservoir coverage must match their
reference partition. `PreparedExpansionSurface` supplies each current area;
electrical length/section never imply an exposed surface. Ambient total
conductance is `h*A_exposed(T)`. A reservoir binding explicitly declares a
uniform areal law, so its native reference `G_ref` becomes
`G_ref*A_reservoir(T)/A_reservoir,ref`. Original boundary events may change h,
reference G and reservoir temperatures; changing coverage requires a new
explicit partition. This does not infer changing bulk conduction lengths,
mechanical support, pressure work or a contact constitutive law.

`resolve_with_boundary_conductances` retains actual area, source identity,
total G and signed heat samples in `ExpandingNetworkSample`. The scalar-only
resolver refuses attached surfaces, preventing a caller from dropping their
thermal feedback. Native circuit, reciprocal field-motion and final-only
entropy paths pass total-G replacements through the original thermal rate and
entropy owners. These are replacements, never additional boundary watts.
Authored reference geometry, reservoir inputs, event cursor and mass survive
acceptance unchanged. The private driver hook moves bounded coefficient maps
once; ordinary drivers use empty maps and the same original numerical path.

Four additional focused G1/G3/G4 groups in
`tests/circuit_material_network/expansion/boundary.rs` compare mixed hot/cold
reservoir and Joule feedback with independent endpoint equations and entropy,
check zero-strain parity and continuing native boundary events, propagate
face-partition-only changes through cooling/current/magnetic motion, and
exercise source/coverage refusals plus late cancellation/replay. Runtime
execution remains separate from test-target compilation and source review.

## Prescribed flow with nonlinear complete-body vibration

`acoustic_realize::realize_blown_ode` prepares either `PreparedExactZoh` for
exactly quadratic original storage or `PreparedNonlinearZoh` for an original
nonlinear remainder. Both own the complete body and viscothermal duct, reuse
all sample buffers, and apply the prescribed flow only at its original inlet.
The nonlinear route retains quartic modal interactions, full cross-mode
resistance, every accepted Maxwell coordinate and the signed surface map.
The quadratic model is not inferred from effort probes, and no nonlinear
stiffness or material branch is removed to meet a render clock.

The common nonlinear owner closes an endpoint discrete gradient around exact
quadratic propagation. Its original-resistance loss is integrated on the
method's discrete effort, and original storage and external work close that
discrete account. This is not exact continuous nonlinear motion or heat; body
and input-band temporal accuracy remain separate from energy acceptance.

Five focused renderer checks include two nonlinear regressions. A warm mixed
body with a 20-microsecond Maxwell branch renders at 48/96/192 kHz. Independent
full-equation RK4 with a separately written physical pressure observer gives
symmetric normalized L2 differences of 2.023823e-5, 5.535631e-6 and 1.495980e-6
on identical physical windows, over 3.6-fold reduction per refinement. Tests
also distinguish the original quartic force and amplitude dependence, retain
warm ringing under zero new flow, transmit an added pulse, keep an unexcited
body silent and preserve all original accepted source history. The existing
complete-state numerical bounds remain; these checks do not establish
arbitrary-clock accuracy, empirical calibration or real-time throughput.

### Spatial interior multiports in the original duct composition

`cavity_modal::interior::multiport` binds original FE boundary modes to one
coherent `GasState`. Lossless, thermoviscous and conductive-wall constructors
require exact FE density/sound-speed agreement with that gas. Transport and
the solid thermal sample use the existing source adapters; no independently
authored air constants or second gas/solid caloric model is introduced.

`duct_scattering` converts the original real-reference FE power scattering
from `exp(+i omega t)` to the duct owner's `exp(-i omega t)` by entrywise
conjugation. Port order, positive real references, modal pressure/flow work
coordinates, and area normalization are unchanged. It exposes the original
row-major matrix for the existing Redheffer owner; callers retain grouping
and matching of connected modes and references. Positive real power-reference
coordinates are not evanescent traveling-wave amplitudes.

`duct_two_port` admits exactly two explicitly retained coordinates and one
common positive real reference, then returns the original
`fs_duct::compose::TwoPort`. Existing `compose_chain`/star composition owns the
network closure, pole/condition handling and terminal load. No retained FE
mode is silently discarded, no feature schema or surrogate fit is added, and
source/cancellation errors publish no partial scattering block.

Focused tests compare actual FE-section plus original analytic duct composition
against one full-length analytic duct, including complex transmission phase,
input impedance and terminal pressure/flow on P2 meshes 1/2/3; finest scaled
complex error must be below 1e-4. Source dry/humid gas and conductive/isothermal
walls change the exported block while the original wall work closes. Typed
source-medium mismatch, unrequested mode truncation, invalid references and
cancellation refuse. These are composition checks, not acceptance of an actual
tone-hole/bend/pad geometry or proof of sufficient evanescent/trace resolution.

## Species-dependent and reacting shock loads on structural surfaces

`WallPressureImpulse::{from_mixture_euler2d, from_mixture_euler3d,
from_reacting_euler2d, from_reacting_euler3d}` consume the actual accepted
per-face transfers of the original evolving-species fields. The complete
reactive owner's report supplies its accepted physical interval; the middle
advection report and final pressure are not substitutes for integrated wall
impulse. Local chemistry changes the gas state without introducing a separate
boundary flux. The source grid retains the physical face area, ordering and
clock, and the caller supplies the explicit orthonormal structural embedding.

All Cartesian pressure adapters share one accepted-face admission. Reflecting
stationary adiabatic slip faces must transfer zero mass for every species
individually, zero energy and zero tangential momentum. Equal-and-opposite
species counterflow cannot masquerade as a sealed wall. Open and periodic
faces refuse support-pressure interpretation. Exterior pressure is subtracted
once over the same accepted area and interval. Existing signed molecular
traction adapters retain their separate `p - tau_nn` semantics.

`WallTractionImpulse::{from_reacting_euler2d, from_reacting_euler3d}` also consume
the original complete molecular/reaction/wave/reaction/molecular interval. Both molecular halves and
the wave stage contribute their actual accepted face transfers, duration and
clock; the signed normal traction can be tensile. It admits inviscid reacting
intervals as the zero-viscous case. The pressure-only reacting adapter refuses
every interval with a molecular report, including compressive net traction,
so a total `p - tau_nn` load cannot be mislabeled thermodynamic pressure.

The original `ShockLoadSurface` projection then preserves nodal force impulse,
resultant, both tangential pressure-footprint moments and affine P1 virtual
work. Four focused test groups exercise a rotated nonuniform 3-D mixture
footprint, planar mixture and reacting supports, actual exothermic 3-D wave
loads that differ from endpoint-pressure estimates, complete fluid/support
momentum balance, boundary refusal and cancelled-step retention/replay.
The reacting molecular case separately retains tensile and compressive faces,
both molecular halves, exterior pressure, P1 virtual work and the previous
accepted flow/operator reports after cancellation.
The volumetric molecular fixture additionally preserves all three support
components and both nonuniform footprint moments under a rotated structural
embedding. Its actual complete face loads differ from the middle wave loads;
positive net traction still refuses the pressure-only adapter. Original
three-dimensional Maxwell--Stefan, heat and full viscous-stress transport and
finite-rate reactions run before projecting those loads. Exterior pressure is
subtracted once, and the original P1 projection retains affine virtual work.
No new gas/source solver, interpolation or inferred geometry correspondence
is introduced. This remains a one-way load transfer from stationary gas walls;
reciprocal moving-wall work, fluid/structure feedback, fracture and material
shock validation require their own admitted physics.

## Authored interior molecular sources in material-string performances

`material_strings` admits `cavity_molecular MAX_CALORIC_COORDINATES MAX_MEMORY_STATES LAW` after the original regional gas records. The existing inline molecular source grammar evaluates only against the selected cavity `GasRecord`, including attached species, NASA-9 and transport sources. The original FE/air-thermal binder selects the joint native caloric owner before preload; exterior observer absorption never supplies the interior source. Complete joint histories replace the heat-only allowance in the existing setup and state budgets. Missing thermal geometry, stale source state, a separate lumped cavity, prepared propagation and incomplete species data refuse before material I/O. Both mean atmospheres remain frozen; these acoustic histories do not evolve mean heat or species inventories. The authored integration checks cover independent interior/exterior gases, original thermochemical source identity, changed native pressure response, and refusals before material reads.

### Source-owned molecular calorics in spatial pressure/heat interiors

`cavity_modal::interior::molecular::MolecularInterior` owns the original
`RelaxationState` and binds its gas density/sound speed bitwise to the original
FE medium. The material evaluates cp once at the requested frequency; that
exp(-i omega t) response is conjugated once into the FE convention. The joint
pressure/heat owner determines molecular compression and gas thermal loss
from the same cp. Total plane-wave absorption is never substituted as an
independent bulk-viscosity or compression law.

The optional conductive wall binds and retains the original
`ThermalTransportSample`, including its source path identity. Common gas/solid
equilibrium temperature, conductivity direction and chosen heat-capacity
constraint retain the existing source adapters' admission. The native prepare
method returns the original caloric FE factor for simultaneous modal flows,
spatial pressure/temperature and separate positive molecular/gas/solid/viscous
work. Material and FE refusals retain their typed underlying cause.

`duct_scattering` exports the full original real-reference matrix by entrywise
conjugation; modal pressure/flow normalization and inward-flow orientation are
unchanged. `duct_two_port` requires exactly two retained coordinates and one
common positive reference, returning the original `fs_duct::compose::TwoPort`.
No evanescent coordinate is silently removed and no fit or auxiliary damping
model is introduced. These remain power coordinates, not conventional
evanescent traveling-wave amplitudes.

Focused source tests use original evaluated reservoirs, independently compare
uniform molecular compression against their bulk modulus, check the duct
phase/conjugation convention, and distinguish zero and changed relaxation
rates. Original prepared solid sources exercise conductivity/capacity-only
changes, preserved source identity, full input/work closure and mismatched
temperature refusal. This is uniform frozen-mean harmonic physics: thermal
mesh/pressure-basis convergence, original thin viscous wall admission and
solid normal-profile limits remain explicit. No mean thermal/species evolution,
nonlinear excitation closure or time-domain realization is inferred.

### Source-owned resolved thermoviscous features in duct networks

`cavity_modal::interior::thermoviscous::source_medium` maps a coherent original
`GasState` to the mixed velocity/pressure/temperature FE coefficients. The
independent bulk viscosity is required explicitly; it must not duplicate a
molecular reservoir supplied by the source. `ThermoviscousInterior` retains
the immutable gas and optional original `RelaxationState`, validates bitwise
agreement with the FE density, pressure, temperature, cp, viscosity and
conductivity, and admits only its declared positive frequency band.

At each frequency, the source heat capacity is evaluated once and conjugated
from the material exp(-i omega t) convention into the FE exp(+i omega t)
temperature equation. Molecular relaxation and heat conduction then use the
same resolved field. No bulk-loss overlay, plane-wave attenuation substitute,
or independent wall resistance is introduced. Responses expose the original
pressure, velocity, temperature and separate physical work terms. Material,
mixed-FE and duct errors preserve their original typed causes.

`admittance` retains the FE convention and one factorization for all columns.
`duct_scattering` conjugates the full direct loaded power-scattering matrix
into the original duct convention, without changing inward flow or coordinate
normalization. `duct_two_port` requires exactly two retained coordinates and
one positive real reference, returning the existing `TwoPort` for
`compose_chain`/Redheffer composition. P is total normal traction; Q is its
work-conjugate inward flow. Nonuniform traction shapes require a matching dual
trace mapping and are not automatically conventional evanescent traveling
waves. No additional port or mode is silently discarded.

Three native integration tests exercise original dry-air and molecular
sources, empty-reservoir parity, changed-rate liveness, source/band/cancellation
refusal, and positive shear/conduction/caloric work with input closure. An
actual resolved feature joins an original analytic duct and complex load at
80/300/1100 Hz; input impedance and terminal pressure/flow match an independent
closure through the FE admittance within 1e-9 of their declared scales. This
is a composition/algebra check, not a continuum or instrument accuracy bound.
The mixed solver's mesh, geometry, equilibrium, thermal-boundary and harmonic
scope limits remain in force. No automatic Helmholtz/LNS mesh interface,
moving mean inventory or nonlinear time-domain excitation loop is inferred.


## Local coating applicability during accepted string performances

`MultiBowPerformance::with_material_shell_coating_monitor` binds an explicitly
budgeted, source-bound real-displacement observer before playback. It checks
the actual initial body coordinates and every accepted mechanical endpoint in
the shared ordinary/contact-thermal step owner, before pressure observation or
resampling. Cached original MITC maps, retained modes, regions, signed film
faces and authored full-tensor strain limits remain unchanged. Detached internal
plans compare immutable source allocation identity without cloning images or
hashing each sample; numeric results, body coordinates and film maxima reuse
fixed admitted storage. No per-endpoint heap allocation or new dynamics solve
is introduced. Preparation and all planned queries are admitted together.

The shared caller cancellation gate reaches local field work. Nonfinite motion,
missing source extents, exhausted resources and an unresolved authored envelope
refuse output. A refused accepted endpoint poisons continuation and retains its
film/sample/bound/limit diagnostics; already accepted physics is not rolled back.
The borrowed report retains original film identity and per-film peak values.
Later preload changes cannot silently replace the checked initial body state.
Bounds cover the represented linear field at observed instants, including mean
motion retained in q; they do not infer missing prestress/cure strain, bound
intersample peaks, certify continuum/omitted-mode behavior, or predict damage.
Nonlinear bow drive limits cannot be inferred from a displacement scale.

Verification: eight actual native Cargo MITC strain tests passed. Four source
plan and three monitor checks compile the exact production owners with original
numerical/cancellation kernels and explicit source/accepted-state containers.
Four exact CLI parser checks and ten launcher tests also pass. These transport
checks do not claim a full coupled-bank, thermal, BEM or violin WAV execution.

## Molecular acoustic heat returns to the same finite source gas

`MaterialWallDg` and `RegionalMaterialWallDg` expose `step_with_finite_gas_molecular_heat_controlled`. Before the scalar acoustic-work residual trials, the original DG trajectory supplies one independently integrated molecular source heat. The existing `ClosedGasParcel` caloric inversion receives external heat plus this internal transfer, while its reciprocal acoustic parameter work remains in the original coupled transaction. The accepted gas then refreshes the same source-derived acoustic operator and original material and molecular memories. Every residual trial starts from the committed owners; failure or late cancellation publishes none of them.

`FiniteGasDgStep::molecular_heat_j` reports the internal source transfer. `exported_acoustic_dissipation_j` is the unchanged native total decrement minus that source integral; roundoff stays visible, and a molecular transfer exceeding the native decrement beyond the declared energy allowance refuses. The combined gas/wave/wall energy residual counts only exported loss and external heat. Generic callback and transport paths retain their previous explicit-export behavior with zero internal molecular transfer.

This path requires the existing `UniformMolecularDiffusion` refresh owner, preserves uniform or spatial radial histories, and retains the fixed-volume, fixed-composition, well-mixed gas and first-order frozen-medium partition. The complete source integral is limited to 128 DG states. Wall heat conduction, viscosity and boundary absorption remain separate exports. The `material_duct_wall` molecular-heat cases cover continuing heating, regional Kelvin/Maxwell memories, independent accounting, source parity and all-owner rollback. New native test execution is pending the unchanged dependency/environment build blocker.

## Material constitutive graph feedback through native spatial enthalpy

`constitutive_graph::spatial::step_candidate` connects original material graph temperature inputs and watt-valued outputs to the heterogeneous spatial enthalpy owner. `GraphSpatialTemperature` samples a convex set of vertices after each vertex has been inverted through its own caloric chart, so unlike specific-enthalpy datums are never averaged. `GraphSpatialHeatPort` normalizes a nonnegative uniform or nodal P1 shape by its native reference-volume integral, then passes the resulting W/m³ field into the original consistent finite-element source assembly. Signed graph watts are deposited exactly once; each output may be bound only once.

The outer unknowns are joule increments on free vertices in ascending order. Every finite-difference and final endpoint trial starts from the committed graph state and original inner coordinates, solves the existing `MaterialGraphResidual`, and evaluates the native nodal enthalpy residual with the selected original radiation and boundary operators. A fully prescribed thermal field still validates and evaluates the graph and bath reaction without inventing a zero-dimensional Newton solve. Explicit vertex, binding, endpoint and aggregate graph-evaluation budgets bound work. Final admission checks every node, the independent P1 heat integral and the full storage/supply/loss account, separating a closed enclosure's numerical net from external supply.

The result carries thermal enthalpy, one graph-history candidate and their physical report together; input histories are borrowed throughout, and refusal/cancellation publishes neither. Fixed reference geometry/mass and backward Euler remain authoritative. This increment binds temperature, duration and heat; it does not add a graph phase-fraction input or phase-dependent electrical law. `constitutive_graph` spatial tests add an independent unit-tetra electrothermal solution, data-only conductivity/heat-capacity substitutions, heterogeneous latent/sensible charts, committed-history checks, prescribed baths and rollback. These new native cases remain unexecuted after the unchanged pinned runtime compiler was killed and shared disk space exhausted.

## Reciprocal native rod and finite-cell enthalpy transaction

`rod_motion::enthalpy::RodThermoelasticSystem` binds every original concentric segment/layer to a distinct whole solid caloric cell. It advances the native rod at accepted local temperatures, then solves the existing finite-cell enthalpy residual with its original conduction and boundary sources. At the unchanged mechanical endpoint it evaluates the full material potential Ψ and entropy S, returns reversible exchange to those same finite stores, refreshes source-derived E(T), G(T) and integrated expansion, and publishes both owners only after all checks pass.

The exchanged heat is `ΔΨ_frozen_motion - Δ(Ψ + T*S)`. The native signed mechanical integration defect is reported separately and is never converted into heat. The total account includes finite caloric storage and the complete mechanical internal energy. Mass/inertia, pose, momentum, supports, attached links, loads and clock retain their original owners. Invalid source domains, duplicate caloric inventory, nonpositive combined caloric/mechanical capacity, residual failure and cancellation leave both histories unchanged.

The admitted model is a first-order frozen-temperature split with fixed reference geometry, fully solid caloric charts, scalar reference density/expansion, scalar or globally affine E/G and the original small-strain envelope. It does not claim phase evolution, section remeshing or damping deposition. Six focused regression groups cover independent six-component potential derivatives, clamped finite heating and material swaps, an axial motion/heat reference, radial and longitudinal conduction, continuous-ODE refinement and rollback. New native execution remains pending after unchanged pinned asupersync compilation was killed and the shared disk filled.


## Source-owned finite solids in resolved thermoviscous acoustics

`cavity_modal::interior::thermoviscous::wall::MaterialWallRegion` and `ThermoviscousInterior::with_conductive_walls` attach the original evaluated `ThermalTransportSample` to each original full-LNS wall patch. Source conductivity direction and caloric constraint remain decisions of the material preparation owner. The adapter uses the existing `conductive_wall_transport` common-temperature admission, retains each original sample through `wall_sources()` in patch order, and forwards original geometry/depth/backing to the FE owner. Gas temperature, velocity and pressure remain resolved, including optional original molecular cp; no gas thermal loss overlay is added.

The source gas is bound before finite walls. `from_gas`/`from_relaxation` refuse FE objects with pre-attached anonymous numeric conductive walls. Source patch attachment is nonempty and once-only; original FE boundary conflicts and frequency-local solid admission remain typed `ThermoviscousError` causes. Immutable borrowing prevents a prepared response from surviving a source change. Existing band, phasor and duct TwoPort semantics are unchanged.

Two new original-material integration cases pass in the focused native fallback: conductivity/capacity substitutions change the full resolved molecular response while all four live loss channels close input work within 2e-8; source identity is retained; wrong common temperature and anonymous numeric wall attachment refuse. The three retained resolved gas/molecular/duct cases also pass. Complete original material, thermochemistry, thermal, port-Hamiltonian and duct libraries compiled in the existing focused seam, with no numerical replacements. This is not a pinned-toolchain/full-workspace gate or an instrument render qualification.

## Prescribed physical bridge loads on the shared string performer

`PhysicalForceProgram<P>` retains physical ports, initial held newtons,
source-timed `ModalForceEvent` assignments, initialization and the original
force resource limits. `StringForceProgram` remains its string-port
specialization. `MaterialForceProgram` accepts `MaterialForcePort::String`
and `MaterialForcePort::ShellBridge { bridge }` in the same authored order.

`MultiBowStringBank::material_force_columns` obtains each bridge row from
the original shell or assembly image and places its complete structural
prefix at the actual flattened body-component offset. All other entries are
zero direct load. A newton amplitude multiplies that original work-conjugate
displacement/load map, including original moment-per-newton lever arms.
The method refuses missing/ambiguous shell ownership, absent bridge indices,
inconsistent extents, nonfinite selected maps, cancellation, size overflow
and the declared projection budget. Sources using residual-flexibility
condensation refuse bridge ports until complete direct-force columns exist;
string ports preserve their original source admission.

`MultiBowPerformance::with_material_forces` and `with_string_forces` use
one binder and the original `ProjectedForceProgram`, released-equilibrium
solver and performance loop. All future events are validated before preload.
Rate and callback limits match the original mechanical clock; controls apply
before the same joint mechanical or thermal step. Duplicate ports sum
coherently before projection, while assignments to one port at one sample
retain the original stable last-assignment behavior. No source projection or
extra physics stepper enters the callback.

`material_force_program` exposes mixed authored inputs;
`string_force_program` exposes legacy string inputs.
`projected_material_forces` and the retained `projected_string_forces`
view expose the same common compiled controls. Preload energy retains the
original full-network account. Finite pulses are held newton loads with
explicit release samples, not implicit impulses, contact impact models or
block-dependent gains. Pressure, material history, coating applicability
and energy diagnostics retain their existing owners and physical limits.


## Source gas and duct composition for hybrid spatial acoustic regions

`cavity_modal::interior::hybrid::HybridInterior` borrows the original source-owned `ThermoviscousInterior` and original pressure `InteriorMultiport`. It verifies pressure-region rho/c and optional boundary-layer viscosity/conductivity/cp/gamma bitwise against the same source gas. Original resolved solid samples remain accessible through `source()`. The source's explicit frequency band applies to every preparation; original FE trace matching and field/work recovery stay in fs-feec. An optional original molecular source now supplies both the pressure region's complex bulk compression and the resolved region's complex heat capacity; the original material is evaluated once per frequency. A molecular source with an existing H boundary-layer configuration refuses because that independent equilibrium thermal correction is not the joint molecular/thermal model.

`prepare` returns the original reusable hybrid preparation. `duct_scattering` conjugates the complete external matrix exactly once from exp(+i omega t) to the existing duct exp(-i omega t) convention. `duct_two_port` requires exactly two external coordinates and uses the original `TwoPort`; no retained mode is dropped. Reference resistances are explicit positive real work-coordinate choices. The original mode/mesh/interface-placement and insulated thermal-interface qualifications are preserved.

Three source integration cases pass in the focused native fallback. Actual original GasState drives both spatial regions and field/work recovery. Hybrid plus original analytic duct and complex load at 80/300/1100 Hz agrees within 2e-9 normalized error with an independent closure using original H impedance and L admittance. Tests pin exact phasor conjugation, source identity, optional H wall transport coherence, and medium/unsupported-wall/band/reference/coordinate/cancellation refusal. The molecular case additionally verifies empty-reservoir parity, original source identity, positive pressure-region bulk and resolved-region caloric work, relaxation-rate liveness and exact duct conjugation. Independent analytic complex-K guide propagation composed with the actual original LNS load agrees to 5.3443575870985403e-8 relative error at the declared 300 Hz fixture, within its 1e-4 discretization gate. This is a source/propagation composition check, not a full narrow-gap or instrument qualification.


## Original molecular source in both hybrid subdomains

The source hybrid preparation now calls the retained `RelaxationState::caloric_response` exactly once per requested frequency. It conjugates cp and K from the material exp(-i omega t) convention into the common FE exp(+i omega t) convention and forwards them to the original hybrid preparation. Original `MaterialError` causes remain typed in `HybridSourceError`. Complete external scattering is still conjugated exactly once when exported to the duct convention. Finite solid samples in the resolved source owner remain unchanged and accessible; K is never overlaid on the resolved heat equation.

The preceding calorically perfect path remains available. Source molecular preparation requires a pressure region without the old H boundary-layer overlay; finite conductive solids remain available in the full resolved region. No mean gas evolution, extra gas loss, loss fitting, time-domain material realization, or full conforming interface is inferred.

The new molecular source case and both retained source-hybrid cases pass on the final implementation. Final focused total across both crates: 25 distinct selected cases, 14 new and 11 retained, all PASS. The coupled final run reused coherent previously verified original dependency rmeta/rlib artifacts with the same Rust 1.100.0-nightly compiler after a redundant material rebuild hit a memory limit. It compiled the current FE and coupled source implementations directly; no production source, dependency or numerical kernel was replaced. Full workspace, pinned-toolchain, DSR/RCH and instrument-render gates remain unexecuted.

## Shared heat network on retained spatial induction histories

`PreparedHeatNetwork::step_spatial_induction_candidate` and
`step_spatial_nonlinear_induction_candidate` place both retained 3D A–phi histories
inside the original complete material heat assembly. Every endpoint trial resolves
thermal-contact resistance, optional integrated layer conduction and enclosure
emissivity/radiosity alongside conductivity and magnetic calorics. Original field
and enthalpy equations, whole-volume receivers and source pins remain the owners;
reports retain the full thermal/interface/radiation samples. Both native/free-energy
and physical internal-energy accounts exclude closed-enclosure numerical net from
external supply. Physical Joule heat and reversible magnetic exchange remain
separate from numerical magnetic dissipation. Graph electrical bindings and
E/J/B/H-dependent thermal boundary laws refuse without their work/state owners.
Fixed geometry, fully solid magnetic receivers and original magnetic source limits
remain. Three new integration groups compare independent two-temperature linear
and saturating references, contact/emissivity substitutions and continuation, then
check field applicability, source-domain exits and cancellation rollback. New
modules pass formatting and whitespace checks; native execution remains unverified
after the unchanged runtime dependency exceeded the authoring worker's memory limit.

### Generated stationary gas on every authored DG wind source

The shared `DgSourceInput` admits `gas_steady CELLS INJECTIONS mole COUNT ID...`
after a `dg_bore` header. The existing inlet, cell, injection and solver records
follow the original cylinder/cone geometry. The same specific input builder now
serves reed and DG sources; it derives cell volumes and lateral wall areas from
the authored bore and invokes the original `SteadyAxialGas` thermal/species
solver. No transport kernel, caloric law or acoustic integrator is duplicated.
The accepted solution's exact `AxialGasField` drives native DG propagation,
viscothermal memory and source inlet feedback. Computed exterior reaction and
microphone observation retain their independently declared common-pressure gas.

Thermal cells and solve controls remain separate from DG mesh/order/wall
quadrature. `DgSourceInput::steady_gas()` exposes the original problem, controls
and solution before preparation. Source, phase, balance, mean-Mach and aggregate
work limits remain active; the existing `CancelGate` reaches geometry and the
original thermal solve. Side openings/chimneys and finite terminal passages or
closures refuse because their steady mass flow and thermal volumes are absent
from this axial model. The field remains stationary throughout playback;
transient gas evolution, condensation and mean-flow acoustics are not added.

Focused native tests compare generated and prescribed fields/records/DG
operators, independent cell geometry, actual flue pressure and callback
partitioning, thermal/species interventions, malformed/source/work refusals and
late thermal cancellation. The bundled flue example reaches computed baffled
radiation, the original microphone and explicit CLI decimation/PCM output.

## Nonlinear spatial materials in continuing circuit and finite-heat solves

`circuit_enthalpy::conduction::circuit::response::ConductionResponseCircuit`
binds the original E/J-dependent conductor laws to the native spatial circuit
owner. Every frozen constitutive trial uses the same accepted charge/flux,
original Kirchhoff descriptor and actual midpoint field. The shared material
response owner evaluates those actual E/J values; its converged candidate
alone advances the electrical history. Initial conditions close the original
material operating point without altering authored charge or flux.

`ConductionThermalNetwork::step_response_circuit_candidate` closes endpoint
temperature/phase, native circuit storage and source response through the
existing finite-enthalpy residual. Spatial E.J and explicitly routed other
resistor heat enter once. Every discarded constitutive trial consumes the
common circuit budget; source/field defects remain numerical residuals. The
candidate publishes thermal state, circuit history and material snapshots
together, with unchanged accepted inputs on source, work or cancellation refusal.

The conduction-enthalpy regressions use independent two-layer thermal balances
and midpoint RC/RLC equations, covering material-only source changes, both
polarities, continuing charge/flux, local heating and refusal/retry. Bounds are
fixed geometry/topology, scalar memoryless conduction and fixed ideal graph
R/L/C. `step_response_contact_circuit_candidate` also retains the original
prepared contact owner inside every nonlinear bulk trial. Actual endpoint
temperatures and phase fractions refresh the same ordered interface sources;
native surface voltage jumps and exact interface Joule heat feed the original
two receiving thermal cells with the authored heat fraction. Interface heat,
bulk E.J and graph-resistor heat enter the same finite-enthalpy residual once.
Independent RLC/two-cell calorimetry regressions cover bulk and interface data
changes, both polarities, accepted history, source-domain and missing-owner
refusal, common work exhaustion and late cancellation/replay. Omitting a
required prepared contact owner refuses; fixed-context bulk-only stepping
cannot silently freeze a contact or discard its heat.

`PreparedHeatNetwork::step_conduction_response_contact_circuit_candidate`
passes the complete original material heat assembly into that same continuing
circuit residual. Thermal-interface resistance, source-integrated layer fluxes
and enclosure radiosity refresh at every shared endpoint, with their original
samples retained in the generic candidate report. Closed-enclosure radiation
is internal exchange; its numerical net receives the existing combined-energy
check. Original E/J-dependent thermal boundary response or applicability refuses
instead of freezing the initial field slice. The focused circuit/heat regression
compares continuing charge/flux, two-sided electrical-contact heat, thermal
contact transfer and radiative exchange against independent energy and RLC
equations, including data-only emissivity and thermal-interface substitutions.

This adds no displacement current, moving material, hysteresis, nonlinear
storage or reversible field-caloric coupling. The first-order coupled thermal
scheme and bounded Picard closure retain their refinement/convergence limits.
### Authored finite-gas heat on the continuing DG aperture performance

`AperturePerformance::with_finite_gas_heat` attaches one original
`ClosedGasParcel` and a finite signed-power program before playback. The first
assignment is at mechanical sample zero; later assignments are strictly
increasing inside the original half-open window. Each held watt value supplies
power times the actual mechanical step to the existing
`ApertureDg::step_with_finite_gas_heat_controlled` transaction. The original
source-supported caloric inversion pays acoustic parameter work and publishes
the same accepted inventory to the continuing aperture, main DG bore and finite
terminal passages. Original wall, diffusion and terminal histories remain in
their owners. The actual full finite-gas step report remains accessible.

Admission requires the same complete initial frozen-gas source, original uniform
circular geometry and full main-gas volume, the original loss/transfer owner,
and explicit caloric/coupling controls. Heat assignments share the existing
performance storage/work budget with pressure, force and boundary programs.
Cold attachment and every physical trial honor cancellation. Failed physical
samples publish neither gas/source/load state nor control cursors and retain
the existing poisoned-renderer contract.


This API retains a fixed-volume closed inventory. It refuses independent
upstream/branch gases and computed exterior mean-pressure jumps; explicit mass
and enthalpy transport remains required for those boundaries.

### Native finite-gas heat programs in plate-valve performances

Native plate-valve files accept `finite_gas`, `gas_energy`, `gas_feedback`,
`gas_heat` and signed `heat SAMPLE POWER_W` records after the complete DG input.
They retain original source packs and derive inventory volume from the native
geometry. WAV metadata records the initial inventory, original physical clock,
power assignments and solve controls before type erasure. This is closed,
fixed-volume gas heating with the existing internal pressure observer. It does
not supply mean mass/enthalpy transport, phase change, wall heating, an
independent upstream/branch inventory, or a computed exterior mean-pressure
jump. Such inputs refuse instead of changing their physical boundary.



### Heterogeneous specimens in one finite moisture atmosphere

`FiniteMoistureAtmosphere::advance_shared` solves one endpoint vapor balance
against all original `MoistureSlab` diffusion owners. Each specimen retains its
own material, geometry, source sorption, optional finite surface-transfer law
and transport controls. Every trial restarts all accepted fields and evaluates
all exposed faces at the same RH. A single gas endpoint receives the aggregate
water transfer; simultaneous uptake and drying do not duplicate the reservoir
or acquire sequential specimen-order bias. Stable IDs define deterministic
assembly order, and the accepted clocks, temperature and total-cell budget must
agree before any native solve. The original single-specimen `advance` uses this
same implementation.

Acceptance checks the represented vapor, every noncondensing species, total
mass and the sum of original solid-water changes. A failed later specimen,
source/domain exit, root failure or cancellation publishes no partial gas or
solid state. Returned slabs remain the original source-bound inputs to swelling,
wet inertia and structural stiffness consumers. Gas composition, pressure,
transport and sound speed come from the original thermochemical owner.

The physical scope remains isothermal, well mixed, fixed gas volume and dilute
bound-water diffusion on fixed dry skeletons. Gas internal-energy change is
reported separately and is not sorption heat or a complete thermostat duty.
There is no condensation, gas momentum or automatic downstream mechanical
publication. G1/G3/G4 regressions cover an independent heterogeneous one-cell
balance over continuing intervals, opposite-sign uptake/drying, canonical input
permutation, a companion material changing an unchanged specimen's native plate
response, and whole-candidate refusal/cancellation.

### Finite-atmosphere sorption calorimetry

The `moisture-sorption-calorics` consumer composes original finite-atmosphere
slab or Cartesian diffusion with `MoistureSorptionEnergy` from the SAME solid
card and exact fixed source context. `advance_with_sorption_energy` and
`advance_cartesian_with_sorption_energy` run the existing conservative endpoint
and then admit every original cell's caloric interval before publishing either
the transport candidate or its energy. Diffusion, finite surface exchange,
the shared RH solve and the species inventory are unchanged physical owners.

At held T, a rigid solid with negligible bound-phase pV stores
`delta_U_bound = m_dry * integral (h_vapor - q_st(M)) dM`. The original retained
water-vapor species supplies the full NASA-9 caloric receipt and molar mass,
even when an active composition has depleted that species. Its datum never
comes from a synthetic liquid-water or latent-heat constant. The fixed-volume
gas+solid system's thermostat duty is
`Q_in = sum m_dry * integral (R_water*T - q_st(M)) dM`.
Positive Q enters the system; negative Q must be removed to maintain T.
Native gas internal-energy change plus bound storage independently audits this
datum-free duty against the caller's absolute joule allowance. Represented
water-balance and gas-energy rounding errors are retained and can refuse the
candidate; arbitrary reference energies never enlarge that allowance.

Local binding/storage increments and mean external thermostat power are
reported. Nonlinear loading-dependent energy changes include internal moisture
redistribution even when net vapor uptake is zero. These per-cell storage
increments are NOT unique spatial heat-deposition sources: internal water-energy
fluxes and conduction require their own closure. Source integration is bounded
by a caller-selected total crossed-segment budget, with no allocation inside
the material chart and cancellation between cells. Only the native candidate
and O(cells) final caloric report are retained.

This rung remains isothermal with fixed dry skeletons and finite ideal vapor;
it does not evolve temperature, add hysteresis or free liquid, or establish
isotherm/caloric thermodynamic consistency or a global entropy inequality.
Focused tests cover native uptake/desorption and the RT correction, sealed
redistribution, actual six-face Cartesian fields, late source/work refusal
and cancellation. Native Rust execution is a separate verification claim.


### Shared material for immersed flow, additional inertia and sphere films

`material_immersed::RetainedFluidSpheres` binds one original `MaterialScaling`
and its `FluidTransport` source to the existing D3Q19 moving-sphere assembly.
The same current transport snapshot determines the bulk BGK relaxation,
reference-density unit conversion and missing-film dynamic viscosity. Film
cutoff/regularization are explicit geometry/calibration inputs; there is no
independent viscosity knob. Each `AdditionalMaterialSphere` retains its
original `ResolvedBody`, requiring a matching uniform spherical material
region. Source density supplies mass and inertia without material-name cases.

`rebind_fluid` admits a new original transport snapshot on that continuing
trajectory and updates bulk relaxation and all normal/tangential film
resistances together. Every accepted population, body pose and momentum,
additional inertia, geometry, physical conversion and step count is preserved.
The original dx/dt and speed bound remain; any changed reference density
refuses because it would reinterpret conserved mass and momentum. Source,
scaling or cancellation failure leaves the original binding and laws intact.
The native collision owner preserves its family and higher-order relaxation;
film geometry, regularization and work limits likewise remain unchanged.

Held SI forces/couples become native impulses using the same dx/dt and density
scales. Actual midpoint gaps drive the original cutoff-subtracted normal film
and optional reciprocal tangent/spin modes. One native solve publishes all
body poses, momenta and fluid populations together. The accepted-step counter
advances only after native acceptance; cancellation/refusal leaves physical
state unchanged. Native receipts retain separate resolved-fluid, film, wall,
external-load and rotational-discretization accounts. Checked conversions
expose their mass, inertia, impulse, angular impulse and energy in SI.

The model deliberately retains the native INTERIOR fluid. Sourced carrier
inertia is ADDITIONAL to it; this is not a displaced-fluid rigid-solid closure,
and subtracting a guessed interior mass is not supported. This lane is
incompressible Newtonian flow with transport held within each step, declared
spherical geometry and stationary
or periodic domain faces. It adds no finite heating, phase transition, contact,
moving-wall model, compressible/rarefied film, measured cutoff calibration or
trajectory accuracy certificate. Low-Mach and native geometric admission are
necessary numerical conditions, not proof of those continuum assumptions.

`tests/material_immersed.rs` supplies three focused G1/G3/G4 groups: original
material substitutions through real bulk/film/body evolution; an independent
SI missing-film impulse/work formula and actual population-plus-carrier
momentum; physical load conversion, continuation and cancellation/refusal
rollback; and source/geometry/allocation-domain admission. Fixtures use
synthetic source cards and do not imply experimental validation.

Three additional focused groups exercise live viscosity intervention after
native motion, exact history preservation and no-op continuation, current
bulk relaxation and independent SI film impulse/heat, original population/body
momentum, density/cancellation rollback with deterministic retry, and native
collision-family/higher-rate preservation. The refresh can consume externally
evolved temperature/composition data but adds no such evolution or finite heat
store; time refinement is still required for varying transport.

## Source-bound hysteretic circuit heating

`ThermalHystereticCircuit` recognizes the original material source retained by
`UniformHystereticCore`. Its receiving core caloric cell must have the exact same
card and complete accepted query. Every fitted magnetic coefficient must cover
the entire declared core temperature interval. Each native enthalpy trial
rechecks those pinned claims at the actual proposed core temperature while
starting from the same accepted magnetic memory. The returned native circuit
retains the accepted source query, flux, polarization, current and clock together
with finite thermal state. Raw authored core fits keep their original authority.

Changing only fitted source data therefore changes actual voltage-driven
current, remanence, core pinning heat and winding Joule heat. The original
Kirchhoff/play and thermal residual owners still perform all numerical work.
Pinning heat and winding heat reach their distinct original cells once each;
numerical residual is not heat. Source-constant coefficients introduce no
reversible field-parameter work. Fixed geometry, constant-density solid calorics,
monotone-linkage steps and first-order thermal coupling remain the model limits.
Temperature-dependent magnetic fits, core spatial eddy loss and moving geometry
are not inferred. Two focused G1/G3/G4 regressions compare independent quadratic
voltage/current and heat balances under a data-only fit change, continue magnetic
memory with actual source temperatures, and exercise card/domain/late-cancel
refusal with unchanged accepted state and deterministic retry.


### Live tensor-solid source continuation

`TensorSolidBody::from_material_sources` binds each Tet4/Tet10 element to its
original complete tensor material query, reference density and material frame.
`at_query_points` refreshes the same pinned claims and rebuilds all equilibrium,
Maxwell and relaxation operators on the original quadrature, support map, basis
and consistent reference inertia. A density change on unchanged geometry refuses
rather than changing the mass of a closed specimen.

`TensorSolidHistory::rebind_materials` preserves displacement, momentum, clock,
physical Maxwell spring strains and accumulated force work/heat. The material
owner recovers spring strain from the old factor and expresses it in the new
factor; copying energy coordinates or resetting a branch is forbidden. Changed
rank-deficient factors refuse at either endpoint: an old nullspace prevents
recovery, and a new nullspace would erase strain needed by later continuation.
Identical factors can retain memory when their relaxation time changes. All six
independent factor directions must pass reconstruction admission, including for
a zero-memory state. The native exact propagator then continues the same history.

The signed storage change at fixed configuration and physical spring strains is
separately accumulated external parameter work and must fit an explicit source
allowance. It is never deposited as material heat. This is prescribed endpoint
coefficient evolution at fixed geometry, not closed thermoelastic feedback,
thermal expansion, phase change, or a guarantee of source support along the path.
No continuous-trajectory or response-band certificate is added. Complete source,
operator and history candidates remain unpublished on every refusal/cancellation.
Three G1/G3/G4 regressions compare six-component factor transfer, live tetrahedral
continuation against independent physical-strain ODE integration, signed work,
native heat, source/mass/domain boundaries and rollback. These are manufactured
constitutive references, not experimental validation.

## Source-retaining material walls on complete DG loads

`DgPortLoad` can retain an original uniform or regional material-wall source
alongside its one native duct, propagator and accepted state. Source-bound
constructors and the authored DG input path preserve these queries, selected
claims, thicknesses, restraints and loss-property roles instead of discarding
them when extracting a duct. Coefficient-only loads do not acquire invented
material sources.

`rebind_wall_temperature_controlled` resolves the original card at one supplied
temperature through the existing pinned-source resolver. Native transfer retains
wall displacement, momentum and admitted Kelvin/Maxwell histories, reports
signed material work, and rebuilds the complete mean-port relation. Original
branches, finite passages, rational terminal poles and accepted clocks remain
with their existing owners. The source, state, propagator and coupling update
publish together only after all admissions and cancellation checks.

The held-flow trial can integrate the original solid-wall dissipation using its
actual solved inlet, terminal and branch-junction inputs. It retains the main
bore's own dissipation allowance before complete-load aggregation; a branch
cannot conceal excessive deposition into a main-wall body. Continuous
trajectories require their own matching source integral and refuse this held
input route. A native heat integral is a physical source report, not permission
to credit a thermal body without a complete energy transaction.

### Multiple finite gas compartments sharing one spatial material wall

`solid_enthalpy::gas::network` generalizes the original cell/gas
backward-Euler residual to independently identified finite compartments.
All existing single-gas entry points delegate to this same residual owner.
Each compartment retains its original material catalog, mass, fixed volume,
composition and optional original implicit reaction mechanism. Stable
compartment IDs define the gas unknown order; globally unique patch IDs and
aggregate exposed-area admission prevent duplicated wall coverage.

Every trial evaluates all gas endpoints, local solid material sources,
source-bound gas/solid layers and the optional prepared material or complete
heat network at the same physical endpoint. Gas patch heat is debited from
the shared wall once, with separate signed transfers and energy residuals
for each compartment. Chemical formation energy remains inside native gas
storage; no independent chemical heat source is added. Original stable
caloric increments are retained instead of subtracting large absolute
formation-energy references. Every gas and wall publishes in one detached
candidate after local and total energy admission.

This is homogeneous, closed, fixed-volume gas chemistry and fixed thermal
geometry. Compartments exchange heat through the solid, not species or
momentum. Returned structural fields retain the existing sequential
thermoelastic interpretation; pressure traction, moving boundaries,
reciprocal thermoelastic heat, flame fronts and a whole-system entropy
theorem are not claimed. Numerical refinement remains necessary for
backward-Euler time accuracy. Source, geometry, work, numerical and
cancellation refusals publish no partial state.

### Source-bound tensor relaxation with finite thermal feedback

`tensor_material::solid::thermal::relaxation` couples original source-bound
tetrahedral tensor dynamics to the original finite enthalpy network.
Every endpoint trial resolves the original material cards at the proposed
caloric temperatures, rebuilds the original material operators on unchanged
geometry/basis/inertia, preserves every accepted Maxwell memory coordinate,
and advances the existing exact-ZOH tensor/PHS trajectory. Independently
integrated element dashpot losses enter their actual caloric receivers once.
The existing common cell-driver, `BackwardEulerStep`, and Newton owner retain
the local thermal residuals, conduction, boundaries, trial budgets, and
fresh final energy admission.

The admitted thermodynamic model has temperature-independent stored elastic
tensors and reference density. Only original positive Maxwell relaxation
times may vary with temperature. Stored coefficients must be validity-wide
scalar constants; equal samples or flat tabulations do not qualify.
Every relaxation curve must continuously cover the complete admitted interval
with positive interior knots. Source validity cannot freeze unowned dynamic
coordinates. The original frequency applicability remains authoritative; this
adapter adds no trajectory-wide frequency qualification.

Mechanical and caloric bindings require the same original card and full query,
unchanged selected density claim, actual borrowed caloric law, and matching
original element mass and volume. Each whole element has one unique thermal
receiver. The caloric chart is fully solid with constant density. Other finite
thermal cells and independent source/conductance callbacks participate in the
same residual; the callback must not duplicate mechanical heat or mutate live
histories.

`RelaxationThermalCandidate` retains one complete proposed thermal/source/
mechanical endpoint. `RelaxingSolid` owns the rebuilt body and the entire native
history; its borrowed history cannot accidentally use another body's factors.
Displacement, momentum, constitutive memory, clock, accepted step count,
cumulative force work, prior parameter work, material heat, and initial energy
survive continuation. A relaxation-only update performs zero parameter work.
Failed source admission, nonlinear iteration, energy checks, or cancellation
publishes none of the candidate.

Thermal rates use backward Euler and endpoint relaxation rates are held during
the original exact mechanical interval. The coupled method is first order in
general; nonlinear tolerances do not certify time accuracy or uniqueness.
Original event cuts and accepted left-limit semantics apply. Geometry,
supports, basis, and reference mass remain fixed. Thermal expansion,
temperature-dependent elastic entropy, phase change, finite deformation, and
basis adaptation are outside this bounded dissipative-kinetics model.

Focused regressions in
`tests/tensor_material_source/thermal_relaxation.rs` compare against independent
physical-strain integration and thermal balance, exercise shared conductive
feedback under companion-material substitution, and preserve accepted state
through unsupported laws, domain exit, context mismatch, and cancellation.

## Finite material-wall heat during coupled aperture performance

`ApertureDg::step_with_finite_wall_heat_controlled` advances the original
moving mechanics and complete downstream load once at the initial material
state. The same held-flow trial can integrate physical Kelvin/Maxwell wall heat
using its actual inlet, terminal and finite-branch junction inputs. The main
bore's own native loss bounds that deposition; unrelated branch, gas, terminal,
mechanical or numerical loss cannot be relabeled as solid heat.

One original enthalpy chart then closes `body storage change + wall parameter
work = external heat + native wall heat`. Every endpoint work trial starts from
the same marched acoustic state and the original body. It refreshes the retained
uniform or regional wall source and the complete mean-port relation without
resetting displacement, momentum, material memories, branches, finite passages,
rational terminal poles or accepted clocks. The moving source receives that
same final load impedance.

The acoustic frame retains its original native dissipation and signed parameter
work. The enclosing body account exports native dissipation minus the heat
already deposited in the solid, measures complete final represented storage,
and checks body, source and combined energy residuals before publishing any
participant. Native terminal exchange remains explicit; no rational storage
function is invented. Refusal and cancellation preserve every accepted owner.

`AperturePerformance::with_finite_wall_heat` owns the original chart, material
card and finite body and schedules signed external power on the physical sample
clock. Heat controls share existing control/work budgets and advance with the
accepted observation. The accepted body and last complete energy account remain
inspectable. Diagnostic observation failure retains the renderer's existing
poisoning semantics.

This is a first-order frozen-material partition for one complete main-bore wall
body, one temperature and caloric chart, an all-solid phase, fixed geometry and
constant density. Regional thicknesses and restraint/source distinctions remain
explicit. Side-branch walls keep their independent materials. Physical time-step
convergence remains separate from exponential/source integration budgets.
Expansion, phase transfer, gas/body conduction, an independently heated upstream
tract and simultaneous finite-gas/wall playback require their own joint
transactions and are not supplied by this route.


### Conservative gas plumbing through one reacting spatial material wall

`solid_enthalpy::gas::network::exchange` composes the original snapshot
nozzle transport with the existing joint gas/reaction/spatial-wall residual.
All nozzles sample the accepted gas network before inventories change.
Species and upstream enthalpy transfer once; every chemical/thermal trial
starts from the identical reconstructed inventories. Finite adiabatic
reservoirs may omit wall patches, while the complete network requires at
least one patch. The complete heat-network entry retains original material
contacts, integrated layers and closed surface radiation.

The enclosing energy account combines native datum-cancelled transport and
joint reaction/heat balances across their identical intermediate state.
Transferred enthalpy and chemical formation energy are never credited again
as heat. Final total gas mass is independently measured. All gas inventories,
thermal history and source-resolved structural fields return together;
source, flow, withdrawal, temperature, numerical, work or cancellation
refusal publishes no partial state.

This is first-order transfer -> joint reaction/heat splitting over one
interval, not simultaneous implicit flow. Accepted wall/chemical pressure
changes affect the next interval's rates; time refinement and event splitting
remain required. Gases remain homogeneous and fixed-volume. Wall geometry
and sequential thermoelastic interpretation retain their existing limits.
Jet thrust, spatial flames, moving boundaries, pressure traction and
reciprocal thermoelastic heat are not claimed.


### Closed affine tensor thermal feedback

`TensorSolidHistory::advance_affine_thermal_candidate` reuses the original
source-operator rebind, native tensor trajectory and shared finite-cell
thermal residual. Caloric exchange equals dashpot loss minus native
coefficient free-energy work minus delta(T*S), allowing signed reversible
cooling. The physical balance uses U=psi+T*S; native mechanical-potential
accounting retains psi and a separate cumulative
`thermal_coefficient_work_j()` account. Thermal coefficient changes never
increment external `parameter_work_j()`, and dissipation enters the thermal
receiver once.

Original card identity, selected claims, complete query context, frame,
exact caloric law, whole-element mass/volume, complete Maxwell memory and
common clock remain required. Geometry, density, supports, eigenstrain and
basis remain fixed; phase change and varying singular factors refuse. The
endpoint-coefficient split is first order and requires physical timestep
refinement. Rejected or cancelled candidates leave accepted mechanical and
thermal histories unchanged. Five focused regressions cover continuing
memory and physical U balance, lossless cooling, full anisotropic rotated
entropy, orthotropic common-scale admission and non-affine/late refusal.

`TensorSolidHistory::advance_material_thermal_candidate` extends the same
driver through `TensorThermalLaw`, including independently affine Young/shear,
Young/Poisson and orthotropic engineering constants with nonlinear stiffness.
One oriented source derivative template
per element supplies exact snapshots at the accepted and trial temperatures;
frozen and evolved states use the same trial snapshot. The shared material
potential supplies mechanical heat capacity as well as psi, S and U. A
whole-source-interval bound must establish positive total caloric plus
mechanical capacity over the original strain domain, even without expansion.
The caloric curve is unchanged, and coefficient work and delta(T*S) still
enter its reversible exchange once. The strict affine entry point retains its
original admission. Focused regressions cover continuing nonlinear material
memory and heat, a held volumetric finite-heat endpoint, and whole-domain
capacity refusal with unchanged accepted inputs.

## Complete material heat on inertial plastic solids

`ThermoplasticDynamics::advance_heat_network_candidate` binds the existing
complete material heat network to the continuing native P1/J2 inertial owner.
Bulk conductivity, integrated thermal layers, source-backed thermal interface
resistance and enclosure emissivity see the same finite caloric cells as local
plasticity and reciprocal thermal expansion. Structural cells retain their
original caloric object, material card and complete fixed source context; a
separately valid thermal query with different pressure or processing coordinates
refuses before an advance.

The original backward-Euler consistent-mass solve retains displacement,
full nodal velocity, plastic/hardening history and the accepted clock. Physical
yield/overstress dissipation and reversible heat enter one thermal evolution.
Elastic/hardening storage, kinetic storage, constitutive numerical loss and
inertial numerical loss remain separate; neither numerical loss is heat.
The original `SolidHeatAudit` removes integrated internal radiation roundoff
from external supply and checks both corrected energy residuals against the
same joule tolerance. Additional finite thermal bodies stay on that shared
timeline. This entry point refuses electrical networks without their work owner.

Continuation accepts returned thermal and mechanical states together and uses
the returned observer's `second` member for independent observer history.
Cancellation, source failure, observer failure and energy rejection publish no
partial state. Geometry, thermal areas/view factors and reference mass remain
fixed; the original fully-solid, small-strain/rotation source ranges apply.
Thermal substeps do not refine mechanical loading or inertia, so physical
interval refinement remains required. No moving contact geometry or melting
claim is added.

G1/G3 checks in
`tests/solid_plastic_material_transport/heat/solid/dynamics.rs` cover the
consistent reference-mass kinetic account, continuing plastic motion, one-time
physical heating, closed total energy and data-only emissivity/interface swaps.
The shared solid test also exercises source-context disagreement and
cancellation for the dynamic entry point. These are synthetic numerical
fixtures, not measured constitutive validation.

## Complete heat networks on reciprocal elastic trajectories

`ThermoelasticSolid::advance_heat_network_candidate` and
`DynamicThermoelasticSolid::advance_heat_network_candidate` compose the existing
prepared material heat observer with the original reciprocal tetrahedral solver.
Bulk conductivity, source-integrated layers, thermal interface resistance and
enclosure emissivity use one prepared heat network; actual equilibrium or
continuing Newmark motion retains the original reversible-heat iteration,
material potentials, nodal work, inertia and finite caloric history.
The shared thermoplastic `SolidHeatAudit` gate excludes integrated closed
radiation roundoff from external supply and checks both corrected residuals
against the original energy tolerance. No thermal or mechanical state is stepped
again for this admission.

Every structural cell must retain its exact caloric object, immutable material
card, material-table binding and complete physical query context in the heat
network. Electrical bindings refuse without their work owner. Candidates publish
the original motion, heat and caller-observer histories together; continuation
retains the returned observer's `second` member. Fixed reference geometry and
supports, the existing small-strain sub-solidus constitutive intervals, and the
explicit thermal timestep/refinement requirements remain unchanged. This does
not supply deformation-dependent contact area/view factors or phase mechanics.

`tests/solid_enthalpy/thermoelastic/heat.rs` defines independent two-cell
thermal/Newmark references, original conductivity/contact/emissivity swaps,
source-integrated conduction, continuing loaded motion, and source/domain/budget/
late-cancellation refusal with deterministic retry. Fixture values are numerical
test data, not measured material validation.

## Reciprocal finite gas and material-wall heat on the continuing DG aperture

The joint finite gas/wall transaction couples one original closed gas parcel
and one original finite main-wall enthalpy body through the existing
MaterialConductionLayers interface. Gas conductivity comes from the resolved
gas internal-energy midpoint. Solid conductivity comes from the original body
card at its enthalpy midpoint, with explicit source conditions, property
selection, area and effective conduction distances. Temperature is reconstructed
from each owner's conserved energy; arithmetic temperature averaging is not
substituted.

The original mechanics and complete acoustic load march once at their accepted
initial state. Native solid heat is integrated once from that interval's solved
ports. Every interface trial returns to the same marched acoustic candidate and
original gas/body inventories. Existing gas-work and wall-work solves refresh
the native endpoint in the declared gas-then-wall order. Their sum is the
represented storage intervention. One signed interface heat is credited to the
gas and debited from the body, and the midpoint conductance law must close its
own residual.

The enclosing balance includes gas internal energy, body enthalpy, mechanical
and acoustic storage, original terminal exchange and external inputs. Deposited
native wall heat is subtracted once from exported dissipation. Interface heat
and both native parameter works cancel internally. The default entry point
keeps gas and terminal dissipation exported; inferred numerical or mechanical
loss is never invented as material heating.

`step_with_finite_gas_wall_viscous_molecular_heat_controlled` opts the original
scalar, paired or modal aperture into the same joint transaction with physical
gas heating. The original solved-port trial supplies the separate positive
bore molecular, bore viscosity and optional solid Kelvin/Maxwell integrals;
main-inventory terminal passages supply their native resolved shear work.
Neither mechanics nor acoustic history is marched again during caloric,
material-work or midpoint-conduction trials. Separate explicit gas integral
budgets cover complete held-flow states of at most 128 native coordinates.

Gas caloric inversion receives its physical loss once, while the finite wall
receives its own selected solid loss. Gas source heat remains distinct from
external heat and reciprocal interface conduction in the retained gas report.
The sum of bore gas and solid deposits is checked against the original main
bore decrement; terminal-channel dissipation cannot subsidize that allowance.
The enclosing exported loss removes the gas and solid deposits exactly once.
Independent inventory, numerical, aperture and unselected solid losses keep
their original ownership. Absence of a gas source needs no integral budget and
preserves the default joint trajectory.

A global material-trial limit supplements the existing individual caloric,
work, interface and energy budgets. Final refusal or cancellation preserves all
accepted inventories, mechanics, material histories, terminal poles and clocks.
Focused regressions cover actual solved bore-flow integrals, independent
terminal shear integration, both caloric balances, absent-source identity,
separate source-budget refusals and final-cancellation rollback with exact retry.
These regressions have source-level checks; native execution remains pending.

The gas subreport describes the gas-first intermediate endpoint, while the
coupled frame describes the final gas-plus-wall endpoint of the same acoustic
march. Their storage must not be added together. This route retains fixed gas
mass/composition/volume, a complete fixed-density all-solid wall body and the
existing first-order frozen acoustic partition. Independent branch inventories,
computed exterior mean-pressure feedback and an upstream tract retain the
closed-inventory refusals. Expansion, phase transfer, resolved thermal profiles,
advection and convection require their own declared owners and are not implied
by the series-conduction model.

## Paired gas and wall heat during aperture performance

AperturePerformance::with_finite_gas_wall_heat retains the original gas parcel,
wall enthalpy chart, immutable material card and material-conduction declaration.
Scheduled signed gas and wall powers share the original physical sample clock.
Each paired event consumes two existing auxiliary control and compilation slots.
The native joint transaction publishes both inventories, mechanics and acoustic
histories together. Heat-control cursors publish with successfully observed
samples; observation failure keeps the existing poisoned-callback semantics.

The accepted gas, wall state, original material/interface declaration and last
complete energy report remain inspectable. Separate gas heating, wall heating
and finite reservoir transport programs cannot own the same evolving inventory
simultaneously. The joint route inherits its fixed-inventory, fixed-geometry,
all-solid scope and the original first-order acoustic partition.

The explicit `with_finite_gas_wall_viscous_molecular_heat` opt-in reuses the
same paired program, physical clock, cancellation gate and original DG owner.
Before any sample, a disposable zero-flow trial admits the original source
and explicit molecular/viscous heat budgets without advancing accepted state.
Every rendered frame then uses the joint source-heat transaction at its
actual solved bore and terminal flows. Bore molecular and viscous heat,
terminal passage viscous heat and optional solid-wall dissipation enter
their declared inventories once, with the original midpoint conduction and
parameter-work accounts. No extra heat-control events or compilation slots
are consumed. The default paired program retains its previous behavior.
The selected budget remains inspectable, and existing observed-sample cursor
publication and poisoned-callback semantics are unchanged.

The focused playback regressions compare direct native joint steps, gas and
wall caloric states, accepted acoustic histories and output across block
sizes, and check early budget refusal and cancellation without cursor drift.
Their sources passed Rust parsing and scoped whitespace checks; native
compilation and execution remain pending.

Focused G0/G1/G3/G4 regression cases in
tests/dg_gas_rebind/finite_performance/wall/gas.rs check independent gas/body
caloric changes, signed midpoint conduction, native parameter work, complete
energy closure, continuing acoustic histories, scheduled powers, block
partition invariance and transactional refusal/cancellation. These cases are
implementation evidence only when executed; they do not supply experimental
validation or establish omitted transport, phase or geometric feedback.

## Authored material walls in the original plate-valve DG performance

The existing plate-valve file accepts optional dg_wall declarations immediately
after its DG geometry/source records. It consumes the shared rigid, uniform,
principal-aligned orthotropic, physical-region and main/chimney assembly grammar.
The original DgSourceInput binds that layout before the actual inlet impedance
and moving plate are constructed. Gas, solid inertia/elasticity, Kelvin/Maxwell
memories, ports, terminal reaction and retained clocks remain in the same native
DG owner. No second wall trajectory or fitted sound filter is introduced.

Wall and plate declarations share one immutable artifact cache. The original
wall resolver bounds its declaration counts and total source bytes; existing
plate-section counts bound additional artifacts. Regional conditions, physical
supports, restraints and thicknesses remain separate. The retained layout is
inspectable on PlateValvePerformance and enters the existing CLI wall metadata.
Omission preserves the original path, and explicit non-DG walls refuse before
material loading. Native geometry, band, source and DOF admissions still apply.

Focused G0/G3/G4 cases in performance/file/wall/tests.rs exercise rigid-path
equivalence, uniform/regional source-to-operator effects, continuing pressure
and wall motion, shared artifact loading, and malformed/source/cancellation
refusals. Tests require execution for runtime evidence; synthetic cards do not
provide experimental material or instrument validation.

## Moving source-bound chambers coupled to shared spatial-wall heat

`cavity_phs::piston::material::solid` connects original material-bound pistons
to the continuing `SolidEnthalpy` gas network. Original-snapshot nozzle flow
precedes joint reaction and finite spatial-wall heat, followed by the native
gas/piston midpoint momentum and pressure-work solve over the SAME thermal
interval. Source bodies, inertial mass, velocity, support and port history are
retained. Final volume, species, temperature and pressure drive the next flow
snapshot. The complete heat-network variant retains original material contacts,
conduction layers and closed surface radiation.

Native datum-cancelled stage residuals telescope across original-to-final gas,
finite-wall and kinetic stores, including exterior pressure, applied force and
support work. Transported enthalpy, chemical formation energy and pressure work
are never credited twice. Final gas masses are independently remeasured; all
moving endpoints, wall history and source-resolved structural state return only
after every stage and enclosing conservation/cancellation check succeeds. The
native substage port report is replaced by one whole-interval tick. Open species
streams explicitly require the enclosing assembly balance, not the closed
piston's two-port audit.

This is first-order flow -> joint reaction/heat -> motion composition, requiring
time refinement and event splitting. Thermal contact geometry and the structural
reference mesh stay fixed. The translating piston is a separate rigid source
body; structural fields retain sequential thermal-loading semantics. Deforming
wall pressure traction, reciprocal thermoelastic heat, jet thrust, crank/contact/
friction laws and spatial combustion are not supplied by this composition.

`tests/solid_reactive_gas/moving.rs` defines independent flow, four-store heat and
midpoint momentum/work references; source conductivity and density changes;
species/formation continuity; continuing port clocks; full heat-network use;
and late mechanical refusal/cancellation with deterministic retry. These native
regressions are authored but unrun in this implementation environment. Separate
synthetic arithmetic checks and syntax review are not native execution or
experimental validation.

## Authored finite wall and reciprocal gas-wall thermal performance

The original plate-valve input now accepts a finite_wall solid block after
optional finite_gas records and before observation. Its wall_caloric record
selects original dimensional heat-capacity and density properties under
single-claim selection and the existing T temperature coordinate. The original
source-derived enthalpy chart retains its claims, query receipts and identity.
The explicit initial specific enthalpy must reconstruct the original wall
temperature. All selected main-bore regions must have the same material card
and condition point; mechanical restraints and thicknesses can remain distinct.

wall_inventory supplies mass, volume and surface area. wall_feedback supplies
the signed parameter-work bracket, source/energy tolerances, iterations and
absolute work allowance. wall_dissipation none or native with dimension, byte
and multiply-term budgets selects the original native solid heat treatment.
wall_heat and wall_power records prescribe bounded signed external powers on
the mechanical sample clock; the final assignment remains held.

A simultaneous finite gas inventory requires explicit gas_wall_conduction and
gas_wall_feedback records. The adapter reuses the original finite-gas parser
and complete-port admitted volume. It validates both original power lists
before merging their sample union, retaining both actual held powers and
charging two shared control/work slots per pair. The existing joint runtime
owns gas/wall energy exchange and both native parameter-work corrections.
Neither inventory is attached independently or advanced twice.

The supported model is one fixed-density all-solid wall body, fixed geometry,
and, when selected, a closed fixed-composition fixed-volume gas inventory.
Original cold/step admissions reject incompatible sources, phases, clocks,
topology, loss owners or resource budgets. Failed/cancelled playback retains
the existing accepted-history and poisoned-callback semantics. No melting,
expansion, advection or resolved convection is inferred. The input grammar is
documented in performance/file/thermal_wall.rs, and the existing CLI records
initial inventories, original chart receipts, paired powers and solve budgets.

Focused authored/manual attachment, unequal-temperature exchange, source
refresh, block partition, budget/refusal and cancellation cases live in
performance/file/wall/tests/thermal.rs. Their runtime assertions require actual
execution; syntax parsing alone is not numerical or experimental validation.

## Reciprocal source expansion through native tensor mechanics and finite heat

`TensorSolidHistory::advance_expanding_thermal_candidate` binds the existing
source-integrated isotropic or oriented-orthotropic expansion model to original
tensor-Prony quadrature and the finite enthalpy residual. Nonzero eigenstrain
adds one stationary unit coordinate to native quadratic PHS storage, supplying
thermal stress through the original exact motion/work/loss owner. It is not an
external force port. Native Maxwell rebind preserves viscous strain through
thermal strain increments, including unchanged singular factor directions.

The original source card, selected expansion claims, stress-free reference and
material/mesh frames remain attached across continuation; generic refresh or
older thermal entry points refuse to drop an accepted expansion binding. Exact
zero expansion keeps the original state extent and trajectory. A stationary
coordinate may drift only within the declared floating-point roundoff budget;
its endpoint is never renormalized or its work/loss recomputed to hide error.

Entropy includes both source-derived stiffness derivatives and equilibrium/
Maxwell expansion terms; the shared finite heat driver receives native dashpot
loss plus signed reversible exchange and checks physical K+U, with U=psi+T*S.
Source-affine alpha(T) may accompany admitted stored tensors over the complete
interval. The existing integrated expansion owner admits each pinned coefficient
knot and transforms alpha/alpha' through its accepted material frame. Numerical
derivatives anchor at the current temperature; the stress-free reference may lie
outside the active interval without extending its affine model there. Actual
free strain still uses the original source integral. Initial, frozen and evolved
states use exact-temperature expansion snapshots, with one shared snapshot for
each trial endpoint. Source-dependent relaxation times remain active.

Actual caloric receiver mass must keep a conservative lower bound on
C_eff=m*dh/dT-T*psi_TT positive on every chart segment. The mechanical curvature
includes .5 e:C_eq'':e, alpha:C_inst:alpha, -2 alpha:C_eq':e,
-2 sum alpha:C_j':q_j, -alpha':C_eq:e and -sum alpha':C_j:q_j.
These expansion-slope terms remain active when alpha crosses zero. Constant
singular branches use the exact factor-memory contraction and an observable
2L envelope, preserving unobserved strain directions without a pseudoinverse.
Independent affine Young/shear sources retain their nonlinear bulk curvature
through the same thermal driver.
The complete interval bound uses the original Mandel elastic-strain limit and
a separately enforced 2L envelope for varying Maxwell spring strain. Accepted,
rebound and final thermal states use that same memory admission. A failed
sufficient bound refuses without claiming that every individual state is
unstable. Constant expansion and the zero expansion law retain their original
paths. Unsupported coefficient paths, expansion derivative jumps, Kelvin thermal
strain rates, finite deformation and phase mechanics remain outside this model.

Five native regressions cover two continuing heated Maxwell intervals, rotated
orthotropic expansion, zero-alpha compatibility, source/capacity/reference/
cancellation refusals, and singular-factor memory. They are authored but unrun.
Pinned syntax/whitespace checks and an independent physical-coordinate RK4 plus
scalar enthalpy reference passed. Actual stationary-coordinate propagation and
its runtime tolerance still require native execution; no experimental authority
is inferred from the synthetic reference.

Three additional regressions exercise continuing affine stiffness/expansion
with retained Maxwell memory, a caloric capacity that passes the old constant-C
check but fails the combined bound, and excessive hidden spring strain.
Independent physical-coordinate RK4 plus enthalpy references were executed;
these new Rust tests were not compiled or run in the unavailable native runtime.
Additional variable-alpha regressions compare two continuing finite-heat
intervals against independent physical-coordinate integration and reject an
insufficient complete-domain capacity even at zero current alpha. Shared
projection regressions check dU/dT and singular observable-memory admission.
These are authored numerical regressions; native execution remains unverified.

## Finite source gas coupled to the actual deforming material wall

`solid_enthalpy::thermoelastic::dynamic::cavity` attaches a closed source gas to
the original inertial tetrahedral wall. Native unit-pressure surface assembly
provides both the nodal force footprint b and affine swept volume b dot delta u.
The accepted gas volume is anchored to the original reference displacement,
rather than accumulated from approximate roots. Initial native equilibrium
loads must equal the same gas/exterior pressure loads plus declared fixed
reference nodal loads; later pressures are generated by the coupled state.

Each bounded Newton-Krylov trial starts from the same original gas, native
Newmark checkpoint, finite-cell thermal history and cloned observer. The two
conservative coordinates are gas delta U and p_start delta V, both in joules.
The shared caloric inversion preserves exact gas mass and cancels formation
datums. Native reciprocal thermoelastic iteration retains local source E(T),
expansion, caloric storage, inertia and exported damping. Independent acceptance
checks enforce geometric volume, gas energy, interface heat, native pressure
work and enclosing gas/wall/kinetic energy, beyond the nonlinear residual norm.

Material-pinned gas/solid layer patches heat the actual bound thermal cells.
Their conductivity follows the original gas and local wall state, with no
second conductivity or heat-capacity law. Patch coverage and ambient-convection
conflicts are checked at binding and every native substep. The original
caloric-slope timestep restriction includes film decay. All caller observer
hooks survive both adapters, including internal conduction, reservoirs,
internal transfers, interval sources, events and timestep limits. Gas-film heat
is removed from external solid heat once, so the same joules leave the gas and
enter the wall. Exterior-pressure and fixed-reference-load work stay external.

This is a well-mixed, fixed-composition, fixed-mass gas with affine small-
deformation volume, stationary supports and reference heat-transfer geometry.
The coupled gas endpoint temperature is held during explicit wall substeps;
this thermal quadrature is first order and requires interval refinement.
Fully solid source restrictions remain with the native wall owner. This path
does not add chemistry, nozzle flow, follower loads, finite-strain volume or
moving phase boundaries.

Three native cavity regressions are authored for independent scalar pressure/
heat agreement in both work directions, large formation datums, eight retained
motion steps with stiffness/density/conductivity substitutions, original
reservoir and interval heat, and failed-candidate rollback. An independent
analytic energy reference closes to roundoff; native Rust compilation and test
execution remain unavailable in this implementation environment.


## Complete heat assemblies on native solid trajectories

`observe_heat_assembly` connects any original `HeatAssembly` to the existing
explicit thermal observer timeline. `MaterialHeatObserver` and
`AssemblyHeatObserver` are aliases over one `HeatSourceObserver`, preserving the
original prepared-network sample and errors. The complete source payload stays
available in the accepted first observer; the second observer retains the
caller's independent history. Providers and samples need not implement Clone.

The four `advance_heat_assembly_candidate` entries use this same observer:

| Original solid owner | Retained native behavior |
| --- | --- |
| `ThermoplasticSolid` | Loaded equilibrium, plastic history and physical plastic heat |
| `ThermoplasticDynamics` | Inertial motion, velocity, plastic history and consistent-mass kinetic storage |
| `ThermoelasticSolid` | Reciprocal equilibrium, material potential, reversible heat and nodal work |
| `DynamicThermoelasticSolid` | Newmark displacement, velocity, acceleration, reciprocal heat and the original clock |

Each entry requires the same original material cards, caloric objects and
physical query context as its prepared heat network. The assembly's finite
coolant cells keep their own original caloric laws. Additional convection,
source conduction, interfaces and radiation use the same actual trial
temperatures and accepted thermal substeps as the solid iteration. Full source
samples survive; no second integrator, mechanical solver or thermal state is
created. Source-bound coolant data can therefore change actual deformation
and retained plastic history through the existing physical response.

Native finite-link conductances enter the original timestep restriction once.
For an explicit observer, `HeatAssembly::explicit_source_conductance_bounds_w_k`
must also supply the original closed-radiation and prescribed-reservoir rate
bounds from that exact source sample. The native convection provider supplies
radiation plus h*A and omits finite links from this additional map. An implicit
provider without declared explicit bounds remains usable by implicit drivers
and refuses explicit stepping. These are frozen-coefficient bounds, not a
nonlinear temporal-error estimate.

The original closed-energy account removes the enclosure's numerical net
from external supply. Finite coolant exchange is internal; prescribed
reservoir heat is external. Physical plastic and reversible heat enter once.
Kinetic storage, plastic algorithmic loss and exported dynamic damping retain
their existing distinct accounts. Late source, cancellation, ownership or
energy failure returns no accepted participant. Continue with the returned
thermal/solid state, the same original assembly and the second observer's
history.

Fixed reference heat geometry, small strain/rotation, sub-solidus source
support and each native mechanical model's domain remain required. Thermal
substeps do not refine mechanical loading; physical intervals need independent
refinement. These entries do not add moving thermal geometry, phase-volume
change, capillary energy or a new constitutive law.



## Source tensor contact storage in the original rigid/finite-heat transaction

`friction_heat::third_body::rigid::tensor::TensorLayerMaterial` now admits the
shared `TensorThermalLaw` potential over the entire receiving caloric
interval. The original material orientation, fixed contact geometry and
`sym(tangent tensor normal)/thickness` strain columns project both stored
stiffness and its complete temperature derivative through the same physical
volume. Equilibrium laws admit fixed-Poisson affine Young modulus, independent
affine bulk/shear moduli, independent affine Young/shear and Young/Poisson
sources, and independently affine orthotropic engineering constants. Shared
whole-interval compliance admission supplies the original rotated tensor,
its derivatives and curvature bounds. The same derivative reaches the retained
contact storage, entropy and finite-heat residual. Independent affine Maxwell
cross terms retain their source identities and original branch order.

The scalar and tensor adapters share the existing `step_heated_impl` force and
enthalpy residual. Tensor entropy is evaluated from retained physical equilibrium
and Maxwell spring strains, so changing rank-deficient arm tensors need no
invented strain reconstruction. Every trial uses the original endpoint-temperature
coefficient convention. Nonlinear sources refresh their derivative at both
actual endpoint temperatures; scalar and affine tensor paths retain cached
constant derivatives. `thermal_slopes()` describes the original temperature,
and `thermal_slopes_at_temperature(T)` supplies an explicit refreshed snapshot.
The caloric exchange is independently computed plastic,
Maxwell and Kelvin loss minus coefficient free-energy change minus
`Delta(T*S)`; total physical storage is `U=Psi+T*S`. Temperature-dependent
Kelvin viscosity and plastic kinetics add no fictitious mechanical entropy.
Actual rigid force/couple work, momenta, retained strain and finite-cell energy
are accepted once, with existing cancellation/refusal rollback.

Whole-layer mass and caloric/source identity remain required. The source-bound
owner admits positive total caloric plus mechanical heat capacity over every
caloric segment using the full original three-dimensional strain envelope.
This sufficient bound can refuse a thermally safe pure-shear projection of a
strongly curved bulk law; no nonlinear derivative is silently frozen.
The fixed-kinematics law refuses unsupported stored-tensor paths, thermal expansion,
changing geometry, phase change and unsupported source intervals. This extends
the retained contact material law; it does not assert measured rosin parameters
or instrument sound validation.


## Distributed material heat in the original bowed contact transaction

The existing `coefficients::heat::MaterialHeatIntegral` now observes the
current unpublished contact trajectory through `evaluate_staged`. It reuses
`PreparedSourceLoss` and its original complete source rows, coefficient
signature and resource limits. Accepted-start displacement, velocity and
Maxwell memory are paired with the actual converged held modal input, including
contact/support effort and preload exactly once. The mechanical owner clears
stale trajectory scratch before preparation, and publishing invalidates the
candidate clock. No second mechanical integrator or loss partition is introduced.

`MaterialRelaxationHeatDestination` binds one complete original material group
to a distinct finite receiving body. The fiber adapter derives its volume from
the original annular layer area, supported length and fibers per width strip;
every original axial/bending mode and Maxwell arm of that layer belongs to the
group. The receiving cell must retain the source card, exact caloric owner,
physical volume and mass, fixed reference density and a valid original
absolute-temperature coordinate. Distinct bodies cannot share one body's heat
capacity, repeat a branch loss, or reuse a separate contact film's cell.

For full-strain fibers, each original strain-history update retains its per-arm
heat before summing. `MaterialHistoryHeatChannel` selects these sources by bound
history identity and original channel index, only from the complete unpublished
candidate of the same sample. The aggregate mechanical law and total-tension
admission remain unchanged. Each destination combines its linear bending
integral with all axial history channels of its physical strip/layer. Duplicate
channels and tokens from independently constructed owners refuse. Preparation
and publication invalidate candidate observation. Existing history laws retain
a default whole-law channel; source laws can expose finer physical channels
without another mechanical update. The v4 contact-thermal profile binds these
same routes through `fiber-body STRIP LAYER CELL_ID` before output creation.

The existing contact/enthalpy Newton candidate deposits independently integrated
body heat once. Original conduction carries it to contact films in that same
trial, where the original temperature-dependent friction law supplies the
mechanical response. The combined account removes that same routed quantity
once from external mechanical loss. Finite heat, contact/source histories,
mechanical state and moving material-coordinate cursors are accepted together.
The existing performance wrapper and moving contact path expose these routes.

This rung admits one uniform temperature per supported strip/layer, a fully
solid body with fixed reference density and geometry, and fixed constitutive
coefficients throughout the source-supported caloric interval. Body storage or
relaxation that depends on temperature, thermal expansion, phase change and
finer spatial body-temperature projection require their corresponding physical
owners. The original linear source integral still refuses unsupported external
radiation states and varying-force trajectories.
The new source/fiber/finite-heat fixtures are manufactured verification cases;
they supply no measured rosin parameters or instrument-sound validation.

### Physical duct observation through material assemblies

`RealizeOptions::baffled_duct` connects ordinary and compiled material-reed
assemblies to the existing regional TMM, moving-reed junction, terminal flow
transfer, retarded Rayleigh microphone and causal `ObserverFilter`. Every axial
segment retains its own coherent material gas; a separately bound, pressure-matched
exterior supplies both radiation loading and observation. The accepted source
wave, including its original priming sample, reaches that same physical receiver.
Internal bore pressure and the compact jet proxy are not included in microphone
output. The returned gas is the exterior; `duct_observer` retains interior states,
receiver geometry and both propagation/filter delay accounts. Solid-reduced reed
mass, stiffness, area and loss remain those of the existing material compilation.

`BaffledDuctObserver::material_wall` binds the existing source-resolved
`LocalMaterialWall` through an explicit `MaterialWallLayout::Uniform` or complete
`Segments` assignment to BOTH the reed's characteristic feedback and the
terminal-flow observation. The local density, elasticity, thickness, axial
restraint, anisotropy and optional Kelvin/Maxwell spectrum are evaluated at each
actual radius and frequency; no wall coefficients or pressure samples are fitted.
The observation retains the complete original wall source. Supplying this binding
beside an authored scenario wall refuses ambiguity. Source frequency support and
thin-membrane geometry admission remain active. The wall's material temperature
is independent of the gas states; neither thermal equilibration nor wall radiation
is inferred. G1/G3 tests compare the independent thin-cylinder reduction, change
only elastic or relaxation data, and check source/geometry refusal.

`AcousticMaterialBindings::duct_wall` selects `DuctWallMaterialBindings::Uniform`
or `Segments` and compiles those sources through the common
material-assembly entry point, including assemblies with material-derived reeds.
The caller supplies the original card/query/claim selection, constitutive source
model, thickness, restraint and frequency applicability. The compiler reuses the
existing wall requirements and reduction, checks both endpoint radii of each
segment, and retains the resulting `MaterialWallLayout`. Realization forwards it
to the same physical feedback and microphone path while preserving the requested
regional gases, exterior medium and receiver. Authored wall coefficients or a
second observer wall override refuse; an absent physical observer also refuses.
G1/G3/G4 regressions compare the direct physical path, swap only elastic or
relaxation source values, retain independent gas/material states, and exercise
unsupported sources, geometry and conflicting bindings. These checks add no
thermal equilibration, shell evolution or experimental validation claim.

`RegionalDuct::with_material_wall_regions` attaches one original wall source to
each authored axial segment, independently of the gas partition. Adjacent runs
coalesce only when BOTH gas and complete wall state agree. Other joins use the
existing pressure/volume-flow continuity and local TMM; material properties are
not averaged and no reflection coefficient is fitted. The retained walls reach
analytic and tabulated terminal loads, stable inlet reflection, characteristic
feedback and physical outlet pressure. Constant-pin shortcuts cannot bypass a
bound source. A second external wall override refuses.

This segment binding is cold and fixed-geometry, with the original 256-segment
cap. Continuous gas fields and tone holes require their existing spatial or
mutual-mouth owner and refuse on this segmented-wall path; the uniform-wall
path retains its previous support. It does not model axial shell transmission,
joint mechanics or thermal evolution. Four added G1/G3/G4 tests compare an
independent two-material Maxwell transfer matrix, uniform/subdivision limits,
source-domain and partition refusals, and a one-segment-only material change
through common compilation and actual microphone pressure. Original wall states
remain independent of interior and exterior gases. These new Rust cases are
unexecuted; nine changed Rust files passed pinned parsing and bounded source/API
review, with scoped whitespace checks passing.

This image's default characteristic route admits unbranched cylinder/cone bores,
their local wall law and a flanged terminal. The explicitly selected lumped route
also admits tone holes under the uniform-gas/rigid-wall scope documented above.
Other body/string/cavity members, extra excitations and unsupported discretization
requests refuse. Existing default assembly paths are unchanged.
Any explicit gas binding requires an explicit observer law, including an explicit
lossless/classical choice with its labelled molecular-relaxation omission. ISO or
molecular laws retain their original admission; scenario air alone preserves the
original default selection. Tests compare original-owner pressure and earliest arrival, independent
interior/exterior changes and a solid-material valve substitution. These are
frozen-media numerical integrations, not mean-flow transport, radiating walls or
experimental instrument validation.


## Source material ownership for spatial finite-wall heat

`ThermoviscousInterior::with_tangentially_isotropic_wall_transport` attaches
complete ordered `MaterialTangentialWallRegion` declarations to already sourced
finite walls. Each declaration supplies the original prepared normal path and
a separately selected tangential path. Admission requires the exact original
normal path/sample, one material-card identity, identical source-state endpoints
and the same density and specific-heat source lines. This preserves the original
caloric constraint and refuses silent cross-material mixtures. In-plane isotropy
and zero normal/tangent coupling are explicit constitutive declarations; they
are never inferred from a scalar normal measurement.

The normal and tangential source samples remain observable. The existing
`prepare_time_with_wall_diffusion` consumes the attached material through the
complete spatial FE wall owner, including all original molecular reservoirs and
structural fields. Normal-only harmonic source responses refuse an attached
spatial wall rather than dropping its material response. `ThermoviscousDg::new`
checks the prepared owner's retained tangential coefficients against those
source samples, so a same-FE normal-only or invented-material preparation cannot
enter the source-bound duct connection. Material replacement requires a new
source assembly; no existing accepted state or prepared coefficients mutate.

This is fixed-equilibrium planar solid heat spreading, with the original common
temperature and physical thickness. It does not add mean-temperature evolution,
phase changes, curved extrusions or a measured instrument-realism claim.

### Distributed native liquid/gas condensation

`liquid_gas::evaporation::multiple::advance_liquid_surfaces_with_gas` connects
several existing native liquid interfaces to the original material finite-gas
chamber. Every surface samples the same gas snapshot and shares its old vapor
inventory. Unique liquid cells, measured adjacent exposed faces, one volatile
identity and one constant liquid/gas datum are required. Canonical coordinate
order selects each actual native mass before caloric reconstruction. Native
samples are subsequently matched by coordinates, independently of sparse-tile
iteration order.

One batched native mass/enthalpy transaction refreshes material state; one native
flow step advances the physical clock. There is no per-patch flow advance or
full-grid copy. The original liquid caloric chart, Stefan-film law, species gas
owner and native VOF transport remain authoritative. Whole total/volatile mass
and source-relative energy accounts include actual native storage and subsequent
flow defects. Absolute local defects also share one interval budget, so opposing
phase streams cannot hide numerical error. Both candidate owners publish only
after source, transaction, flow, balance and final cancellation admission.

`advance_liquid_surfaces_with_gas_and_heat` adds the same shared material
chamber's simultaneous reciprocal sensible heat before joint phase exchange.
Every native patch supplies an explicit coefficient and admitted liquid heat
branch. Exact postheat liquid states enter the phase owner, and the native
storage account retains sensible heat carried by mass leaving after heating.
Combined source heat/phase defects are counted once before native storage and
flow defects. Phase, heat and combined reports remain distinct, while one native
transaction and one flow step store the full original-to-final enthalpy change.
Zero heat preserves the existing phase-only endpoints. Three focused regression
sources exercise mixed heat/phase counterflow, zero-heat parity and heat-only
exchange; the mixed case resolves the heated-mass carry above its energy budget.

These phase/flow and heat/phase/flow splits are first order and require complete
timestep refinement. Existing single-patch entry points retain their behavior.
The routes retain positive-interface capacity,
fixed connected gas volume and original electrical/surface-history admission;
they add no dry-wall nucleation, gas recoil,
sealed-vessel displacement or overflow/exhaustion remap. Four regression sources
cover cold/warm counterflow, independent whole-system mass/energy accounting,
input-order invariance, one-clock continuation, singleton equivalence and
transactional refusal. These sources do not replace native Rust execution.

### A finite thermal wall with a zero-capable liquid film

`finite_gas_thermal::wall_film` couples the original fully-wetting material
boundary to the existing uniform-body thermal and finite-gas caloric owners.
A cold dry wall can acquire condensate, warm as it receives condensation heat,
and later lose its entire film. The wall and liquid retain separate masses,
caloric charts and energy storage; a dry surface supplies no liquid capacity.

Nested bounded scalar feedback closes the wall endpoint and reciprocal
sensible heat. Each trial heats the original gas, applies the original phase
boundary at the trial wall temperature, and checks convection against the
actual final gas. The original body solver consumes the same phase heat debit
and retains its conductivity, radiation enclosure and Biot admission. Film and
wall temperatures must agree at acceptance. Radiation remains an external
source, and vapor-carried sensible/latent energy is not applied a second time.

The returned uniform-body callback carries candidate wall, gas and film states
for one full requested interval, even when its frozen-rate film dries earlier.
Source-addition and boundary-port rounding enter the signed and absolute
energy accounts. Whole-system and volatile mass changes are independently
checked. Failure or cancellation publishes no partial physical state.

This is a first-order, fixed-area, fixed-gas-volume, fully-wetting planar
closure requiring complete-timestep refinement. Original uniform-wall
admission remains necessary; a small convective Biot number alone does not
justify localized rapid phase heating. Nucleation barriers, drainage, liquid
motion, displacement work and a coupled entropy theorem are outside this
boundary. Four focused regression sources check new-film heat capacity,
complete dry-out against an independent coupled caloric solution, the original
sensible-only limit, and atomic refusal/cancellation.

### Continuing wall-film formation and drying

`WallFilmThermalEnvironment::advance_window` repeats the same original
wall/film/gas transaction on one bounded physical clock. The supplied body and
gas heat are totals for the entire requested interval, applied at uniform rates.
Cumulative time endpoints apportion those sources, including the last partial
interval. Source-apportionment roundoff is independently charged to the common
energy account. Rejected phase-step proposals consume work but contribute no
physical heat, mass, event or elapsed time. Explicit step-limit refusals retry
from the same accepted state; nonlinear/source refusals propagate unchanged.

One attempt ceiling and whole-window absolute energy/mass allowances apply.
Accepted shorter intervals resume the original nominal time-grid endpoint.
Actual final wall storage, total and volatile mass, and carrier inventories are
checked independently. Compact samples retain physical temperatures, actual
film mass, onset and within-window dry-out events. The caller publishes wall,
gas and film together only after the complete window succeeds.

This subdivision satisfies admission limits; it is not a temporal-error
estimator. The maximum physical step must still be refined. Five regression
sources cover uniform source totals, donor-driven retries and global budgets,
a continuous dry/wet/dry/wet lifecycle, a live four-rung full-coupling time ladder,
and refusal after an already accepted internal interval. Native test execution
remains necessary before claiming measured convergence or validated physics.

### Condensate storage in finite gas/body heat exchange

`finite_gas_thermal::condensation` consumes the original source-bound
`fs-material::phase::condensation::WetGasParcel` and the existing uniform-body
enthalpy callback. `heat_candidate` solves the actual combined gas/liquid
internal-energy increment against supplied joules, using the same bounded
scalar feedback owner as other finite thermal couplings. Initial nonequilibrium
water may change phase and temperature even at zero external heat. Film
evaporation is limited by its actual inventory; every carrier species remains
fixed. The material owner aligns the NASA/IAPWS caloric reference once, and its
reference-invariant sensible-plus-latent increment drives the thermal root.

`CondensingGasThermalEnvironment` closes one backward-Euler reciprocal transfer
between that wet parcel and the original finite body. Radiation belongs to its
separate enclosure. The convection callback sees the accepted candidate's
actual vapor depletion, source-selected transport and temperature. Reports
retain the phase/species account, the body trajectory, source-law residual and
measured body + gas + liquid storage residual. The enclosing joule limit must
pass independently; a latent-heat diagnostic is never reapplied as an external
source. No accepted input is mutated on success, error or cancellation: the
caller publishes both candidate physical states together.

This is homogeneous near-atmospheric liquid-water equilibrium within the
material owner's explicit temperature/pressure bounds. The retained liquid
shares gas temperature and does not displace its prescribed volume. Wall-film
kinetics, droplets, ice, pressure work, gas flow, nucleation and a coupled
entropy/passivity theorem are outside this lane. The original body conductivity,
source-domain, radiation and Biot restrictions still apply. Focused G0/G3/G4
tests cover dew-point crossing, finite-film reheating, continuing finite-body
cooling, source-aware convection, mass/energy conservation and refusal rollback.
## Finite solution activity, bound-water transport and native swelling

`moisture::solution::FiniteMoistureSolution` retains the shared
`fs-material::osmosis::IdealDiluteSolution` inventory. Its independent solute
amounts remain unchanged while water uptake concentrates the solution and
desorption dilutes it. Explicit solvent molar mass and partial molar volume
define the incompressible water mass scale; no material name supplies a density.

The equal-pressure, isothermal equilibrium-vapor bridge maps aqueous water
activity to the original `DryBasisSorption` relative-humidity coordinate, under
the same pure-water reference pressure and ideal-vapor approximation.
Every scalar trial advances the original `MoistureSlab` from its accepted field
using the same endpoint activity at all exposed faces. Original diffusivity,
surface-exchange, source-domain and solver controls remain active. A bounded
safeguarded secant solve closes solution solvent change plus bound-water change
over an explicitly supplied source-supported volume bracket. Only bracket
residuals are retained between trials, avoiding two additional full slab and
material candidates. No inventory is clipped, and dissolved solute is not lost
or diluted by changing concentration independently of its amount.

The returned reservoir, native moisture field and original source receipts
publish together after both incremental and direct water-inventory checks.
`MoistureSolid::advance_finite_solution` resolves original source-integrated
hygroscopic strain, compatible deformation, stress, wet inertia and current
density at that same endpoint. A late mechanical/source refusal preserves both
accepted water inventories. Existing native plate resolution consumes the same
field and refreshes its actual elastic sources and wet mass.

G1 manufactured one-cell uptake/desorption endpoints independently check water
balance and equilibrium activity; G3 tests cover finite dilution, unchanged
solute amounts, continuation, sealed boundaries, equilibrium and salinity's
effect on native plate/solid responses. G4 tests cover inadequate volume and
iteration bounds, temperature mismatch, cancellation and a late native elastic
source refusal. These authored cases do not establish material calibration.

This bridge assumes nonvolatile solutes, solute exclusion from the bound-water
solid, and negligible vapor-path storage. Its original vapor surface resistance
does not become a liquid-immersion permeability. It does not add a resolved vapor
space, ion transport, sorption heat, stress-driven transport or pressure-driven
liquid filtration. Reported osmotic pressure remains a thermodynamic diagnostic
and is not silently substituted for a mechanical surface traction. The dilute
thermodynamic and original fixed-reference moisture limits remain explicit.
Mechanical response is one way: elastic work does not feed solvent chemical
potential, and the route does not claim closed chemo-mechanical energy balance.

`FiniteMoistureSolution::with_solubility` explicitly retains a
`FiniteSolubilitySolution` instead of assuming every component remains
dissolved. Every scalar volume trial repartitions the same total stocks, then
uses only dissolved particles for endpoint activity and the native sorption
query. A surviving solid buffers its solute concentration; crossing saturation
can create a solid, and dilution can exhaust it. Accepted candidates preserve
the complete phase owner for continuation, with independent dissolved-plus-solid
component balances before publication. The shared controller carries this
behavior through single and multiple slabs, Cartesian fields and their existing
mechanical consumers without inventing a second moisture integrator.

Focused G1/G3/G4 cases check manufactured buffered uptake/desorption,
precipitation and solid exhaustion during real native diffusion, retained
crystals on continuation, the resulting native wet-plate mass difference and
failed-solve immutability. This opt-in route has the original material law's
independent neutral-solute, instantaneous-equilibrium scope; it does not claim
electrolyte solubility, crystal-water transport or dissolution heat feedback.

### Finite-rate crystal stocks in the native solution/moisture bridge

`FiniteMoistureSolution::with_kinetic_solubility` retains the original
`KineticSolubilitySolution`, its conductances and growth-surface policies.
At each candidate endpoint solvent volume, the material first preserves both
actual phases and then advances them for the proposed physical interval.
The resulting dissolved-particle activity drives the existing source queries,
native diffusion fields and shared scalar water balance. Every rejected root
trial restarts the original phases; no trial can cumulatively dissolve the bath.

The existing pure-solvent shortcut is disabled for kinetic inventories:
zero initially dissolved material with a positive crystal stock can acquire
solute and change activity during the interval. Sealed specimens still permit
internal solution phase evolution at fixed solvent volume. Nonvolatile solute
and crystals remain in the reservoir; no crystal-bound water or chemical
solvent production is introduced. Accepted endpoints retain identical original
kinetic sources and component totals, with an explicit dissolved-plus-solid
balance before publication. Cancellation maps to the native moisture refusal
and leaves both the reservoir and every specimen unchanged.

Single slabs, Cartesian fields and shared collections reuse the same original
endpoint controller. Freezing reservoir volume at the endpoint is a first-order
approximation to the continuously coupled phase/water trajectory. Root
convergence enforces water conservation, not temporal accuracy. Existing
adaptive wrappers compare moisture, solvent mass and activity, without separate
per-component phase tolerances; phase composition may require a finer physical
interval. Temperature remains prescribed, with no phase heat feedback or solute
transport into specimens. This does not add a second diffusion integrator.

Five focused regressions exercise a manufactured initially-undissolved
dissolution/uptake endpoint, a closed transfer surface, exact sealed-volume
kinetics and continuation, seeded/unseeded growth, original Cartesian/shared
field updates with opposing water transfers and stable reordering, and rejected
roots/work/cancellation preserving original inventories.

### Simultaneous specimens sharing one finite solution

`FiniteMoistureSolution::advance_shared` and `advance_shared_cartesian` couple
independently sourced slabs or Cartesian specimens to one actual aqueous
inventory. They reuse the existing `MoistureParticipant`,
`CartesianMoistureParticipant` and `SharedMoistureSolid` views, re-exported by
`moisture::solution`. Both single-specimen APIs now delegate to the same shared
scalar controller. No sequential bath debit or separate activity root is
introduced per specimen.

Every root trial starts all original fields at one common endpoint activity.
Each specimen retains its original geometry, dry mass, source diffusivity,
sorption and optional surface resistance. Stable distinct IDs set evaluation
and output order; accepted clocks and solution/source temperatures must agree.
The aggregate cell budget includes sealed specimens, whose internal diffusion
still advances. One bounded safeguarded solvent-volume solve retains the
original ideal or nonideal thermodynamic model and every nonvolatile solute
amount. Rejected trials retain only scalar residuals between root iterations.

Compensated local water changes preserve signed uptake and desorption even
when their sum is nearly zero. A balanced net bath exchange therefore accepts
the independently evolved specimen fields instead of treating them as locally
equilibrated. Pure solvent keeps activity exactly one and debits the sum of
actual specimen exchanges once. Final incremental and independently accumulated
stored-water balances must both meet the same fixed coupled tolerance; that
allowance is not multiplied by participant count. Solution and all native field
outputs publish atomically after source, inventory, work and cancellation
checks. These simultaneous methods advance one caller-specified physical
interval; their shared adaptive entry points below complete a whole window.

Focused authored tests cover manufactured simultaneous uptake/desorption,
nonzero opposing exchanges at zero net bath depletion, permutation/replay,
nonideal six-face Cartesian exchange with a sealed diffusing companion,
pure-solvent and single-participant limits, and late source, cell/root budget,
duplicate-ID, clock and cancellation refusals. These are numerical/source
checks, not measured-material calibration or a closed thermal/chemical-work
balance; the original equal-pressure, negligible-vapor-storage bridge remains.


## Consistent contact mobility on finite solids

`finite_solid::impulse::FiniteSolidPointMobility` exposes the complete row-major
point/direction response `B M_ff^-1 B^T` on the original consistent reference
mass. `point_mobility_matrix_per_kg` retains cross-point and cross-direction
entries, component-specific stationary supports and each bounded native mass
solve report. The row cap bounds its dense output; the original mass controls
bound each solve. Cached point traces and one response column at a time avoid
retaining a mesh-sized response per contact. No lumping or symmetry repair is
performed; contact consumers still admit the actual matrix and response.

`point_impulses_candidate` accumulates all material-point impulses through the
same P1 transpose map before one native mass solve, including repeated points.
The original angular/support/work acceptance gates therefore retain all joint
work terms. The single-point API delegates to this shared path. Tests use the
analytic consistent tetrahedral inverse, including negative cross-vertex
mobility, support elimination, combined impulse work and endpoint agreement.

## Rotating material-point loads on native rod motion

`rod_motion::PreparedRodLoads` binds world-held forces at original material
points, world-held nodal forces/couples and material-following free couples.
`RodMotion::step_spatial_candidate` evaluates each nodal resultant at both
actual endpoint frames through the original momentum-kick, centroid-drift,
free-rotation and accounting owner. A rotating force arm contributes `(R r) x F`
at each kick; it is never frozen to its starting world moment. The existing
`step_candidate` fixed-world wrench API uses the same physical stage.

The original kick/drift/kick interval retains six low momentum words per node
beside the candidate's native high words. A first half-impulse therefore reaches
both the COM displacement and the second kick even when it is below the spacing
of the stored momentum. The same retained spin drives swept-angle admission,
the original compensated `fs-time` free-rotation owner and the endpoint kick.
For full inertia, that owner conjugates the increment through the fixed
principal frame without a rounded full-state basis roundtrip. Applied kick work
pairs each impulse with the retained mean velocity through the same mass and
inverse-inertia map. Existing position remainders and eccentric COM offsets are
included before the translated position rounds.

This scratch is initialized after the first material-loss flow and discarded
before the final one. Torsional, bending, axial and transverse loss operators
consume the rounded native state; no old impulse tail can reappear after
dissipation. Accepted state remains f64, supports retain their original masks
and reactions, and the whole candidate keeps its original refusal semantics.
Native regression definitions cover signed force and axial torque impulses,
uniform source-owned rods, and a displacement exposed by large drift
cancellation. This does not change the temporal splitting order, make the
complete timestep allocation-free, or claim exact finite-step energy.

Offsets and reported nodal moments refer to material nodes. The existing native
kick alone transports them to eccentric mass centroids and retains the full
central tensor, support reactions, power and orbital/spin momentum accounting.
Cold preparation reuses the mass-independent compensated `fs-mbd` first-moment
reducer, with at most 16,384 points across the rod. Per-step load evaluation has
constant storage and work per node, independent of attachment count.

Applied work remains the native symmetric kick quadrature. Conservative held
point forces have their usual potential; general free couples do not acquire a
potential or a symplecticity claim. These interval-held loads do not implement
velocity-dependent drag/friction, closed-loop motors, contact or moving supports.
Focused fixtures cover scalar nonlinear torsion and work refinement, rotating
body versus world couples, independent point power/frame covariance, full-tensor
eccentric COM loading, supported moments and atomic refusal/replay.

`RodThermoelasticSystem::step_spatial` carries these same prepared loads through
the original frozen-temperature mechanics, finite-cell enthalpy solve and source
refresh. Both thermal and mechanical histories publish only after all original
material and energy gates pass. Applied force/couple work and the signed native
integration defect retain separate accounts; numerical error is never deposited
as heat. The complete thermoelastic split remains first order. Regressions check
fixed-world parity, independent rotating attachment impulses, complete material
plus caloric storage, and late thermal refusal followed by a valid retry.

## Intrinsic material torsional loss on native rods

`rod_motion::PreparedRodTorsionalLoss` binds the existing
`ConcentricTorsionalSection` shear-viscosity law to every original rod segment.
The enriched loss bundle may add an explicit viscosity property, but its card,
complete physical query, selected original mechanical properties, annular
geometry and shear convention must match the actual model. Transient admission
requires fixed-state coefficients with explicit broadband frequency support.
A temperature/material refresh makes an old loss source refuse; the caller must
prepare the matching new source rather than retain stale viscosity.

`RodMotion::step_torsional_candidate` and
`step_spatial_torsional_candidate` reuse the existing motion owner under fixed
world or prepared rotating-point loads. For relative log rotation phi, the
viscous rate is d(phi_x)/dt. Its exact left/right SO(3) Jacobian pullbacks retain
finite-bend coupling and equal/opposite world moments. Resistance is
`c = sum(eta_s J)/L`; no loss ratio or measured calibration is invented.

Each frozen-pose two-node relaxation uses the original full central inertia,
with supported rotations removed from mobility. An exact exponential impulse
produces nonnegative shear heat and the actual support reaction. Bounded stack
scratch and a palindromic segment sweep avoid a dense damping solve. Symmetric
half-loss/native-step/half-loss composition is second order for the admitted
smooth fixed-state model. Each loss map is dissipative even for stiff viscosity;
elastic wave stability and timestep refinement remain required.

The native receipt reports total and per-segment physical heat separately.
Its numerical defect is delta mechanical energy plus heat minus applied-work
quadrature. Heat is not silently deposited into a reservoir, and this path adds
no finite thermoelastic reciprocity, bending/axial viscosity, contact, fluid-film
loss or temperature evolution. Focused synthetic fixtures cover exact pair heat
and reaction under load, strong viscosity, two-segment analytic damped modes and
heat refinement, finite-bend gradient probes, collective spin, full-tensor
eccentric/frame momentum balance, source mismatch, thermal refresh and atomic
refusal. They do not establish experimental material validity.

## Intrinsic bending loss on native rods

`rod_motion::PreparedRodBendingLoss` binds the existing material-owned
`ConcentricKelvinVoigtSection` to each original rod segment. It retains the
source geometry, card identity, complete query and original mechanical-property
selections, requires explicit broadband fixed-state support, and uses
`C_b/L = sum(eta I)/L` for both bending directions. Missing viscosity is not an
explicit zero. An updated physical source requires a newly admitted loss law.

`RodMotion::step_damped_candidate` and `step_spatial_damped_candidate` accept
`RodAngularLosses` with independently optional bending and torsional laws.
For `phi = Log(R_i^-1 R_j)`, their respective rate components are
`(phi_dot_y, phi_dot_z)` and `phi_dot_x`. They share the exact left/right SO(3)
Jacobian pullback, original central inertia, exponential scalar relaxation and
support-reaction owner. No finite-rotation spin-difference approximation or
second damping implementation is introduced. Bending's two directions and
adjacent segments form a bounded palindromic sweep; simultaneous torsion and
bending compose symmetrically. The original torsional-only API/flow remains.

The native step reports positive bending and torsional heat separately in
original segment order. Its signed numerical defect includes each physical
loss once. Frozen-pose loss is passive for stiff viscosity; second-order
composition still requires elastic-wave stability and temporal refinement.
No dense damping matrix, per-edge heap allocation, inferred axial/shear loss,
bending memory, thermal deposition or thermal bending feedback is supplied.

Focused synthetic regressions cover independently differentiated finite-angle
elastic motion and exact pair heat in both bending planes, strong viscosity,
small-angle analytic damped-motion and heat refinement, explicit-zero and
pure-twist isolation, finite-rotation frame objectivity with eccentric inertia,
simultaneous separate loss accounts, and atomic source/band/step refusals.
They do not establish measured material calibration or spatial convergence.

## Intrinsic axial loss on native rods

`rod_motion::PreparedRodAxialLoss` binds one original
`ConcentricAxialViscosity` per segment to either prepared rod inertia view.
It preserves the original elastic sources, complete query, property selections
and reference geometry, and requires explicit broadband viscosity support.
Its resistance is `C_a/L0 = sum(eta_a A)/L0`; missing or frequency-local loss is
not inferred as a transient law. Temperature/source refresh requires a newly
admitted matching loss bundle.

`step_axially_damped_candidate` and `step_spatial_axially_damped_candidate`
compose axial loss with independently optional `RodAngularLosses`. They use
the native analytic derivative of geodesic-midpoint extension `L0*gamma_x`,
including frame dependence, through original nodal masses and complete central
inertia. Material-node covectors are pulled back to the same eccentric COM/spin
coordinates as native motion. Rigid translation and common spin annihilate the
strain rate, while constrained components receive their actual support impulse.
The shared exponential scalar relaxation supplies positive integrated heat;
bounded stack scratch and symmetric segment/mechanism sweeps avoid a dense
damping solve and retain second-order splitting for the admitted smooth model.

`RodMotionStep::axial_heat` retains total and per-original-segment physical heat,
separate from bending, torsion and signed numerical error. Every physical loss
appears once in the step energy gate. Explicit zero viscosity preserves existing
motion and angular-loss arithmetic. Refused or cancelled trials publish no
state. This mechanical path exports heat without implicit thermal deposition,
and supplies no transverse-shear viscosity, hereditary axial spectrum or
temperature feedback. Elastic-wave stability and temporal refinement still
apply. Focused fixtures cover exact split motion/heat, stiff loss and reactions,
ringdown refinement, finite-pose derivatives, frame objectivity, source matching
and atomic refusal; they do not establish measured material calibration.

## Intrinsic transverse shear loss on native rods

The native fixed-state rod also admits transverse shear through
`PreparedRodTransverseLoss` and `RodMaterialLosses`. The section retains its
original shear correction and explicit per-layer shear viscosities, giving
`C_s = kappa * sum(eta_s A)`. The native geodesic-midpoint derivatives of
`L0 * [gamma_y, gamma_z]` supply both rate and force/couple impulses, including
finite-rotation and eccentric-inertia terms. Axial and transverse loss share
one displacement-flow owner; bending and torsional loss remain independent.
Exact frozen-pose exponential flows compose symmetrically across directions,
segments and mechanisms. They report separate transverse heat, support
reactions and the signed integration defect before publishing a candidate.
Explicit zero viscosity preserves the preexisting trajectory. Sources must
match the original mechanics and supply broadband fixed-state support.
This Kelvin section model does not resolve cross-sectional warping or hereditary
shear relaxation. The fixed-state native step owns no thermal reservoir; the
finite-cell connection below supplies explicit thermal feedback. Elastic-wave and
splitting accuracy still require timestep refinement. Focused G0/G1/G3/G4
tests exercise manufactured two-direction decay and independent heat, stiff
viscosity, a three-clock refinement ladder, zero-law parity, rigid motion,
source mismatch and cancellation. They do not establish measured calibration.

## Finite-cell heating from intrinsic transverse rod shear

`RodThermoelasticSystem::step_material_viscous_with_transverse` and its spatial
counterpart add an independent `RodTransverseViscosity` selection to the existing
axial, bending and torsional choices. The compatible
`step_heat_assembly_with_transverse` and spatial counterpart carry the same
selection through original contact, radiation, coolant and reservoir transport.
Every path uses one native mechanical advance and the same caloric transaction;
the original selection structs and existing entrypoint arguments are retained.

Each segment's actual two-direction shear heat is distributed by its original
`kappa * eta_i * A_i` coefficients into the original whole-layer cells. The
section supplies the shear correction exactly once. Layer heat capacity affects
the resulting temperature, not the dissipation weights. Scaled partitioning
retains representable small deposits at extreme source contrast, and both
segment and total heat closure are gated before publication. The interval's
sources stay frozen during thermal trials. The SAME selected claims must
support the accepted endpoint; subsequent intervals resolve the heated state.

`RodThermoelasticReport::transverse` retains source sections, destinations and
actual deposits, with independent `transverse_generation_j_k = sum(Q/T_end)`
under the existing endpoint thermal quadrature. This internal conversion enters
the finite stores once and is excluded from external thermal supply. Reversible
exchange, other physical losses and numerical defects keep separate accounts.
Source exhaustion, cancellation or a late energy/assembly refusal leaves both
accepted histories unchanged. Focused regressions cover unequal layer capacities,
independent oscillator heat, heated source refresh, simultaneous loss mechanisms,
world/spatial loading, zero-law parity and endpoint-refusal retry. No hereditary
shear law, phase transition or deformation-dependent heat geometry is added.

## Finite-cell heating from intrinsic shaft viscosity

`RodThermoelasticSystem::step_torsional` and `step_spatial_torsional` connect
the original source-bound torsional loss flow to the existing whole-layer
finite-cell enthalpy owner. An explicit `RodTorsionalViscosity` identifies the
shear-viscosity property and one claim-selection policy per original layer.
The existing material extension API preserves every mechanical receipt and the
complete accepted typed query while resolving viscosity from that same card.
Temperature-dependent viscosity is freshly resolved at each accepted step.

Mechanics runs once with accepted-temperature coefficients. Each segment's
physical heat is partitioned by its actual layer eta_s J contributions and added
only to those already-bound whole caloric cells. This source and partition stay
frozen during nonlinear thermal trials. Explicit zero viscosity deposits zero;
no link, unbound store, imposed source work or numerical defect is heated.
The receipt retains the used section sources, per-layer viscosity/polar moment,
native mechanical heat, actual deposited heat and their independently gated
partition residual.

True external thermal input excludes both signed reversible transfer and the
actual intrinsic-loss deposits. Complete mechanical plus caloric storage thus
retains the original native integration defect as its separate numerical
account. Nonnegative `sum(Q_layer/T_end)` is reported separately as endpoint
torsional entropy-production quadrature under the existing first-order split,
not as an exact continuous integral of dissipation divided by temperature.
Reversible, thermoelastic temporal and caloric temporal accounts remain distinct.

After the original endpoint mechanical-source refresh, the SAME selected
viscosity claims are checked at the final actual material state. Another claim
cannot conceal an exhausted interval source domain. Both histories publish only
after that check and all existing energy/caloric gates; failure permits an
unchanged-state retry. This remains the original fixed-geometry, solid-layer
thermoelastic model, without phase change, fluid-film heat or a new time solver.
Focused manufactured fixtures check unequal-layer partition and temperature rise,
changing eta(T) across continuing steps, energy/entropy closure, zero-loss parity,
prepared spatial loads, and late endpoint-source refusal.

## Finite-cell heating from simultaneous rod bending and torsion

`RodThermoelasticSystem::step_angular_viscous` and
`step_spatial_angular_viscous` connect the existing native bending-loss owner to
the same accepted rod and whole-layer caloric history. `RodAngularViscosity`
selects bending and torsional viscosity independently. Each enabled mechanism
resolves its explicit property from the original accepted material state; an
absent selection does not assert a lossless material. Explicit zero retains a
zero physical heat account.

The single native motion step retains both angular-loss flows and their separate
per-segment heat receipts. Bending heat is distributed by the original material
owner's layer `eta_b I` contributions, common to both circular bending axes;
torsion retains its independent `eta_s J` partition. The original fixed cell
mapping, rather than cell mass or a temperature average, determines every heat
destination. Each bending segment and the total deposited heat must match the
native loss within the declared coupling tolerance before publication.

Both deposits enter the original enthalpy residual exactly once, alongside the
signed reversible thermoelastic exchange. Physical bending loss is excluded
from external supply and reported as a separate nonnegative endpoint
`sum(Q_bending/T_end)` entropy-production quadrature. Native integration error,
caloric temporal error and thermoelastic temporal error retain their original
accounts. This is the existing first-order thermal split, not an exact
continuous entropy integral or unconditional time-accuracy claim.

The next accepted interval resolves viscosity from its actual heated state.
The current interval's selected bending and torsional claims must also support
the complete thermal endpoint. Source, thermal, energy or cancellation refusal
leaves both original histories unchanged. This adds no material law, loss
ratio, phase transition, deformation-dependent heat geometry or second clock.

The existing `step_heat_assembly` and `step_spatial_heat_assembly` entrypoints
accept the same independent angular selections while retaining compatibility
with `Option<RodTorsionalViscosity>`. Complete source contact, radiation,
finite-coolant and reservoir transport therefore shares the original motion,
finite-cell energy solve and atomic publication with both intrinsic losses.

## Finite-cell heating from intrinsic axial rod loss

`RodThermoelasticSystem::step_material_viscous` and
`step_spatial_material_viscous` carry axial, bending and torsional viscosity
through the same native rod interval and caloric transaction.
`RodMaterialViscosity` selects the mechanisms independently;
`RodAxialViscosity` identifies the original axial property and one source claim
selection per layer. Existing angular and torsional-only entrypoints retain
their semantics. Complete heat-assembly calls accept either the new selection
or the existing angular/torsional forms.

The original material section supplies both the axial mechanical coefficient
and its ordered layer contributions `eta_i A_i`. Each segment's actual native
heat is partitioned by those same stored contributions, with independently
gated segment and total closure. Reference areas, original material receipts,
selected isotropic parameter bindings and finite-cell destinations are retained.
Cell mass, heat capacity and bending inertia never substitute for the loss
weights. Different finite heat capacities determine the subsequent layer
temperature rises through the original nonlinear enthalpy solve.

Mechanics runs once at accepted temperatures. The axial source and partition
stay frozen through thermal trials; endpoint source refresh preserves native
pose, momentum, mass and clock. The SAME selected viscosity claims must support
that endpoint before publication. The next interval resolves eta at its actual
heated state. Explicit zero viscosity deposits zero. Missing or exhausted
source domains, cancellation, failed heat partition or any late physical gate
leave both accepted histories unchanged.

`RodThermoelasticReport::axial` retains sources, layer areas/viscosities and
actual mechanical/deposited heat. `axial_generation_j_k` separately reports
nonnegative `sum(Q_axial/T_end)`, an endpoint production quadrature under the
existing first-order thermal split. Axial heat enters the finite stores once
and is excluded from true external supply. Reversible exchange, physical
angular heat and signed numerical defects retain their original accounts.
This adds no alternate thermal integrator, hereditary axial law, transverse
shear loss, phase change or deformation-dependent heat-transfer geometry.
Additional-link conductances remain live assembly outputs at every thermal
trial; their cold admission values are never substituted for evolving physics.

## Compensated translation on native rod trajectories

Native rod drift retains one private world-coordinate translation remainder
per material node. The existing error-free `fs_math::eft::two_sum` operations
carry displacement below one coordinate ULP into subsequent intervals, instead
of aborting valid bending/torsion continuation or discarding the same small
motion repeatedly. Eccentric free rotation also retains roundoff in the
centroid-to-material-node transport. Public poses and force evaluations remain
at their rounded f64 coordinates; this compensation does not supply arbitrary
precision geometry or erase the original numerical integration defect.

Stationary supported coordinates do not drift. Original momentum updates,
support reactions, finite/overflow admission, clock representability and energy
gates remain authoritative. Clone/restart, instantaneous impulses and source
refreshes retain the remainder with the same trajectory. A failed or cancelled
candidate cannot consume it. The original reciprocal bending/torsion-strike
continuation and a uniform-translation regression at one eighth of a coordinate
ULP per interval verify physical continuation, accumulation, clone/restart,
retention across an actual velocity impulse, and refusal followed by retry.

## Point impulses on continuing geometrically exact rods

`rod_motion::RodMotion::apply_point_impulses` applies simultaneous terminal or
attachment impulses to the existing centered or eccentric trajectory. Each
`impulse::RodPointImpulse` identifies a material node, a point offset in its
rotating material frame, a world linear impulse and an optional world free
couple. The original source-derived central inertia supplies the angular
response. The point arm is shifted by the actual mass-centroid offset, retaining
translation/rotation coupling, products of inertia and world spin plus orbital
angular momentum.

The operation changes actual nodal momenta without advancing pose, time, elastic
motion, cached forces or accepted step count. Duplicate-node impulses accumulate
before the jump, retaining their frozen-mass work cross terms. Stationary
translations and frame clamps receive balancing linear impulses and pure couples
about the material node while doing no work. The existing refusal of eccentric
free rotors with constrained material-node translation remains authoritative.

`impulse::RodImpulseControls` requires positive finite absolute tolerances for
linear impulse, angular impulse and work. Acceptance checks each node and the
complete batch; opposite unresolved transfers cannot cancel into an accepted
global account. Actual stored linear-momentum and central-spin increments
determine world angular change about the supplied origin. The kinetic increment
uses those actual increments dotted with mean velocities, and the work residual
subtracts impulses before multiplication. An unresolved transfer returns
`RodMotionError::ImpulseBalance` without publishing the candidate. The report
retains these stable changes, caller-applied impulses and support reactions.

Cost is O(nodes + impulses) with O(nodes) scratch, using the existing prelowered
inertia. Cancellation is checked per event and node. The caller supplies the
physical strike, actuator or contact law and any loss partition. Tests compare
an eccentric velocity jump with an independent six-by-six spatial-mass solve,
check clamped/sliding support couples, propagate a reciprocal strike into native
bending/torsion, and reject sub-resolution parallel kicks while retaining the
stable energy of a resolved transverse kick against large baseline momentum.

## Reciprocal isolated impacts between native rods and finite rigid bodies

`rod_motion::RodMotion::apply_rigid_normal_impact` resolves one authored normal
contact between an original rod material attachment and an original finite
`fs_mbd::RigidBodyState`. It derives the rod's point mobility from its actual
mass, central tensor, centroid offset and supported degrees of freedom. The
rigid participant supplies its native full-inertia point mobility. The existing
`fs_contact::impact::solve_normal_impact` selects the restitution impulse; the
original rod and rigid impulse owners apply equal and opposite actions. Both
continuing states are returned together, retaining poses, elastic/material
history, clocks and step counts. Inactive or separating witnesses preserve the
original state bits and require no native impulse event.

The witnesses must have a nonnegative signed gap and share the declared normal
line within the geometry bound. Their actual separation still contributes to
the angular-momentum admission check. Stationary rod supports receive their
original physical force impulses and pure couples, with zero support work.
Actual stored momentum changes, endpoint relative velocity, and pair kinetic
work must satisfy independent absolute bounds. Double-double sums of the
momentum changes retain small collision losses during large common translation;
the relative speed subtracts common translation before adding point rotation.
An unresolved linear impulse cannot be hidden by a resolved angular response.
Failure or cancellation publishes neither participant, including a final pair
balance refusal after both native impulse candidates have been prepared.

The original contact report's dissipation is mechanical restitution loss,
without an inferred heat partition or damage law. This is a frozen-pose isolated
frictionless impact, not contact discovery, a simultaneous contact solver,
persistent contact, penetration repair or an additional rod integrator. Existing
rod support restrictions remain binding. Contact/mobility preparation is
constant work; the original rod impulse path is O(nodes). Focused tests compare
collinear finite masses with the closed-form collision law, offset/full-tensor
response with an independent six-dimensional spatial-mass solve, rotated and
translated systems, support reactions, continuing motion, small loss during
large common translation, and refused/cancelled events with a valid retry.

### Coulomb glancing impacts on the same continuing states

`RodMotion::apply_rigid_frictional_impact` extends this interface through
`rigid_impact::friction`. The original centered/eccentric rod inertia, supported
translation/frame masks, rigid full tensor and original material arms supply
one complete normal-first 3x3 Delassus matrix. The existing
`fs_contact::frictional_multicontact::solve_frictional_contacts` solves restitution
and both Coulomb tangent components together using an explicitly selected
interface coefficient. One accepted world vector acts on the rod and its
opposite on the finite body through their original impulse owners. A nonzero
tangent impulse requires collocated witnesses within the geometry allowance;
their actual separation still participates in angular admission, without moving
either point.

The normal and frictional entries share the same compensated audit of actual
stored momentum and pair kinetic increments, original stationary-support
reactions, and independent endpoint/work gates. Both returned states retain
original clocks, histories and continuing integrators; zero impulses preserve
their bits. Coulomb determinacy, budget, source-state, geometry, cancellation
or late balance refusal publishes neither state. Total impact loss remains
distinct from a thermal/friction partition. Focused tests cover analytic
sticking/sliding, zero-friction parity, independent 6x6 eccentric-inertia
response, world-frame/boost/tangent-basis covariance, additional glancing spin,
stationary support reactions, continuation and transactional retry. This
remains an isolated frozen-pose event, not persistent contact, collision
discovery or a simultaneous-contact island.

## Simultaneous normal contact networks of original spatial plates

`plate_impact::network::PlateImpactNetwork` extends the original spatial plate
impact owners to a graph of independently meshed plates. Each
`PlateNetworkContact` identifies two original owners and one shared world-x/y
point, with a declared restitution coefficient and activation distance.
`PlateImpactSurface` retains each owner's original `PlateDynamics`, complete
consistent mass, support map and explicit reference face height. A positive
edge impulse acts along world +z on its positive endpoint and -z on its negative
endpoint. Mixed load implementations use `L = dyn PlateLoad`; no load, mass,
stiffness, physical damping or trajectory owner is substituted.

Preparation reuses `PlateRigidImpactSet` for each referenced owner's incident
points. The one shared matrix is
`W_ij = sum_p s_pi s_pj B_pi M_p^-1 B_pj^T`, where absent incidences contribute
zero and `M_p` is that owner's original scaled reduced consistent mass. All
cross-interface entries on the same plate remain, even when contacts lie on
different interfaces or one edge initially separates. The cached matrix is
assembled once with compensated accumulation. An event delegates the whole
admissible normal law to `fs-contact::multicontact`; neither independent edge
solves nor a sequential sweep approximates the simultaneous impact.

Actual gaps combine both native interpolated displacements with the declared
reference face separation. Every referenced owner must have a finite matching
checkpoint at exactly the same native physical time, including owners reached
only by inactive edges. Unreferenced states need not share that time and are
not admitted or copied. The original contact law rejects penetration, retains
its complementarity and determinacy results, and rejects globally active
energy gain rather than silently reducing restitution.

An accepted event accumulates each owner's signed transpose request in raw
double-double full nodal coordinates and rounds only the reduced impulse sent
to the original native update. At most one original structural impulse solve
and acceleration restart runs per affected plate. The common original plate
transfer audit reconstructs the actual full momentum change, retains
constrained-row cross inertia, and computes each fixed-row support impulse as
actual momentum minus the unrounded requested impulse. Every free row must
meet its original impulse bound. A complete exactly zero nodal request does
not restart its owner; an all-zero contact event leaves all states bitwise
unchanged and returns no native receipts or balance record.

Before publication, actual endpoint velocities must match the shared contact
prediction, actual relative impulse work must match the original law, and all
owners' unrounded momentum and kinetic increments enter one collective audit.
Only represented transverse Pz and bending Lx/Ly are claimed; complete support
reactions are external. Actual aggregate kinetic change equals contact work
and balances the original mechanical dissipation within declared bounds.
Per-owner reported scalars are rounded only after this common audit. Reports
preserve original owner indexing, including absent receipts for untouched
owners, and original edge indexing for gaps, speeds and law quantities.
Velocity and acceleration publish only after every gate and the final
cancellation poll; positions, clocks, accepted-step counters and native
telemetry remain in place. Refusal, native load failure or cancellation before
publication changes no trajectory.

Admission caps the original owner count (at most 128), global contact count
(at most 64), total full nodal dimension and combined live workspace before
graph allocation. Each original owner retains its own point-location,
full-mass, contact-count, cross-mobility-product and native workspace budgets.
The combined bound includes referenced-owner workspace ceilings, graph arrays,
compensated requests and staged receipts; borrowed model/trajectory storage and
the contact solver's independently bounded scratch are excluded. Work uses
bounded cancellation checkpoints. Preparation assembles signed cache entries
in `sum_p k_p^2` work after original point-response preparation; an event
aggregates each owner's incident requests and performs at most one native
impulse solve for that owner. These are structural complexity statements, not
measured throughput claims.

Focused G0/G1/G3/G4 tests use an independent three-mass centroid-chain solution,
owner and edge permutations, a sequential-sweep falsifier, continuation through
the original damped time owners, nonmatching spatial interpolation and dense
consistent-mass inverses, original support cross inertia, redundant contact
determinacy, equivalence to the existing two-plate owner, coupled initially
separating and inactive rows, global energy-gain refusal, exact no-op behavior,
budget admission, synchronized clocks, heterogeneous load restart/failure,
cancellation and deterministic retry. This is fixed-normal instantaneous
contact. Collision discovery, persistent-contact stepping, rotating normals,
in-plane/drilling momentum, tangential friction, classical impact adjoints,
thermal loss destinations, experimental validation and measured performance
are not established by this extension. A unique velocity witness does not
establish unique local pressure, support load shares or heat partition.

## Native rod-to-spatial-plate simultaneous normal impacts

With `rod-dynamics`, `plate_impact::rods::RodPlateImpact` exchanges momentum
between one continuing original `RodMotion<I>` and one original spatial
`PlateDynamics<L>` trajectory. Each `RodPlatePoint` names an existing material
node and its rotating material offset; the matching prepared plate site uses
the original interpolation and its transpose. Mesh x/y are world coordinates,
the plate normal is world +z, and an explicit reference face height is added to
the original transverse displacement. Both native clocks must agree. The rod
receives +z impulse and the plate its reciprocal -z action.
Original world positions and rotated material offsets remain compensated
through alignment and gap evaluation, so large common origins cannot round
small penetration into an inactive separating witness.
The same native point-kinematics and directional-mobility helpers serve the
existing rod/rigid contact graphs. Their original endpoint identities remain
intact while origin subtraction and material offsets precede final rounding;
both signed normal gaps and tangential misalignment use that shared geometry.

The original cached plate block `B M^-1 B^T` is augmented at each event by the
actual rod's supported cross mobility. Centered and eccentric rods share their
native mass, full central tensor, centroid offsets and stationary support
policy. Shared-node rotational/translation cross terms are retained, while
distinct rod nodes still interact through the original consistent plate mass.
One global `fs-contact` normal law resolves all active rows together, including
initially separating rows. Its nonpenetration, work admission and
unique-impulse/unique-velocity distinction remain authoritative. Redundant
impulses do not establish unique local pressure or a heat partition.

One original rod point-impulse batch and one original plate structural impulse
plus acceleration restart produce private candidates. Actual original point
speeds, each native stored momentum transfer, both support resultants and the
complete event's momentum/work must pass before either participant changes.
Original rod sources, pose, elastic response, translation remainder, clocks and
step counts remain with the trajectory. The plate retains its displacement,
clock, accepted steps and telemetry. Background damping acts through its
original acceleration restart and subsequent native integration; the impact
does not turn damping into an additional impulse or force divided by a timestep.

Compensated native momentum and kinetic increments are combined before final
rounding, retaining small collision work beneath large common translation.
Rod support force/couple reactions and plate constrained-row cross inertia are
included about world zero. All represented linear and angular increments must
close: the transverse plate contributes Pz and bending Lx/Ly, while the rod
retains its full spatial accounting. Failure, cancellation, mismatched clocks,
misaligned points or unresolved stored transfers publish neither candidate;
an all-zero event invokes no native update.

Preparation and event storage have explicit contact, plate mass-product,
original rod-node and combined workspace bounds. The normal law admits at most
64 sites. Preparation reuses the original consistent-mass solve per site;
events add O(k²) rod cross responses, O(k n_plate) interpolation/transpose work,
one O(n_rod + k) native rod batch and one native plate mass solve/restart.
There is no per-contact trajectory clone, surrogate rigid plate, invented
inertia, contact-time search, moving-normal law, persistent-contact stepper,
tangential friction, shell in-plane dynamics or mechanical-loss thermal route.

Five focused physical groups in `tests/plate_impact/rods.rs` compare a full
six-by-six rod spatial-mass inverse and dense original plate mass against
three-site simultaneous impacts; check permitted frame/translation and contact
permutation invariance, restitution endpoints and continuation through both
native integrators; reconstruct stationary rod/plate support reactions; compare
distinct-node and redundant contact motion; preserve small loss beneath common
translation energy; and exercise late native-load/cancellation refusal,
geometry, clocks, budgets, exact rollback, retry and zero-event identity.
The geometry refusal includes a -0.25 m gap under a 2^53 m common origin that
would otherwise round to a false +0.25 m gap before contact admission.
The shared rod/rigid graph regressions retain positive separating gaps and
reject negative penetrating gaps under large translated origins. The combined
original-source native build passes all 104 tests: 32 plate/contact tests,
including all five rod/plate groups, and 72 original rod tests. This selected
native closure is not a full workspace Cargo or DSR verification claim.

## Native sustained normal contact between two spatial plates

`fs_couple::plate_sustained::PlateSustainedPair` advances two continuing native
linear `PlateDynamics` states through one finite contact interval. It connects
the original `fs_contact::sustained::solve` force/displacement natural map to
the actual endpoint gaps produced by the original operator generalized-alpha
integrator with `rho_inf = 1`. Trial contact forces act through the original
point transpose rows, positively on the positive plate and negatively on the
negative plate. Consistent mass, linear stiffness and both Rayleigh damping
terms retain their original native owners. The force response does not use an
impact impulse divided by a timestep, a penalty spring, a virtual contact mass
or a reduced surrogate trajectory.

The lane has one fixed world +z normal, common world xy contact points and
explicit undeformed face heights. Both native states must start at the same
finite time, with every declared gap closed and every relative normal speed
zero within the caller's dimensional tolerances. Initially open contact,
penetration and impact activation are refused. Contact forces are held during
the interval. The actual native endpoint q/v values define the declared
quadratic midpoint gap interpolant; its analytical extrema must remain
nonpenetrating, and a contact with force above the numerical force tolerance
must remain closed with zero endpoint normal speed. Zero-force detachment is
allowed. An open endpoint requires detached dynamics or an independently
located subsequent activation before another call to this lane. This does not
claim continuous-event detection, moving normals, tangential mechanics or
contact-force uniqueness.

Every original mechanical response starts again from the same accepted q/v and
recomputes its initial acceleration under that trial's held contact forces.
Thus a previous active set cannot silently supply the wrong initial force.
Complete external generalized loads are supplied at start, midpoint and end,
including all fixed rows; every retained row must match the original load
callback exactly at its corresponding time. A held-load constructor repeats one
sample. Final accepted accelerations are independently restarted under the end
source and accepted contact force. The source and Rayleigh work receipts are
explicitly midpoint quadratures; no exact continuous forcing or damping
integral is asserted between samples.

Original complete mass and stiffness matrices must be finite and symmetric and
match every retained entry of the native reduced operators. The independent
final audit uses the actual stored q/v changes to check midpoint kinematics and
each free generalized equation. It recovers stationary-support impulses from
the complete rows of M, C and K, including constrained-row cross inertia and
the complete source/contact loads. The current step's internal elastic force
and moment must have zero resultant within tolerance; foundation forces are
not silently treated as internal. Transverse Pz and bending Lx/Ly are the world
momenta represented by this plate model.

The reported balance is native momentum change plus modeled damping impulse,
less external and support impulses. In particular, mass-proportional Rayleigh
damping exchanges momentum with the modeled environment; it is not advertised
as internally momentum-conserving friction. Kinetic and elastic energy changes
are formed directly from actual increments and mean states. Each plate must
separately satisfy its midpoint work identity, damping quadrature must be
nonnegative within tolerance, reciprocal contact work must agree with the
actual gap change and ideal normal-contact work must vanish. Raw double-double
quantities are retained through the combined momentum/work gates. The receipt
separates damping quadrature, source quadrature, contact work and numerical
energy residual. None is automatically deposited as heat, and no numerical
damping parameter below `rho_inf = 1` is admitted by this lane.

Contact count and natural-map workspace/response limits belong to the original
force solver. Point/mass preparation, complete stiffness, native solver
iterations/scratch and trajectory histories are separately bounded here. Trial
states contain one native step receipt rather than a copy of the accepted
history. The original law's mandatory final mechanical replay supplies the
candidate; no extra unbudgeted force response is inserted after convergence.
Both native trajectories publish together only after the path, individual and
combined mechanical checks, physical endpoint restarts, history reservation
and final cancellation check. On refusal the accepted q/v/a, time, step count
and history entries remain unchanged; reservation may change Vec capacity.

Four new focused native integration tests pass: repeated loaded closed contact
against the analytic consistent-mass acceleration and force, including stale
initial acceleration and continuing histories; common-motion Rayleigh decay
against the midpoint factor with approximately fourfold endpoint-error reduction
under timestep halving; an independent dense 13-by-13 supported elastic
midpoint saddle system with both Rayleigh terms and complete reactions; and
mode/source/workspace/history/cancellation refusals, zero-force detachment,
native detached continuation and guaranteed late endpoint-restart rollback.
The focused native runs cover 108 distinct integration tests in total (104
impact/rod tests plus these four sustained tests), compiled against the
original selected source owners. This is scoped native validation, not a full
workspace Cargo test claim.

## Simultaneous normal contact graphs on native rods and finite bodies

`rod_motion::RodMotion::apply_rigid_normal_contact_graph` resolves authored
contacts between one continuing rod and a bounded slice of original finite
`fs_mbd::RigidBodyState` participants. An oriented `graph::RodRigidGraphContact`
can join a rod material attachment to a body, two distinct bodies, or two
distinct rod nodes. These topologies can coexist in one original
`fs_contact::multicontact` solve. A contact cannot join two attachments belonging
to the same rigid participant. Body indices identify the caller's original
participant order; rod offsets are measured from the material node, and body
arms from the original center of mass.

Each signed row retains both physical endpoints. The common graph owner
assembles the complete `J M^-1 J^T` response from original native central
tensors, centroid offsets, attachment arms and admitted stationary supports.
Shared participants supply cross terms with their actual endpoint signs;
unrelated mass blocks contribute nothing. No surrogate bodies, diagonalized
inertias or pairwise impulse sequence replace that response. Shared translation
is subtracted in double-double arithmetic before adding point rotation. Finite
normal directions are normalized with a scaled norm; actual witness gaps and
tangential separation remain authoritative. Penetration is refused.

The original contact owner supplies globally admissible restitution impulses,
passivity, total mechanical loss and the unique-impulse versus unique-velocity
classification. One original rod point-impulse batch and one original resultant
momentum update per affected body prepare the candidates. The same represented
world vector acts oppositely at each edge's two endpoints. Compensated force and
couple sums retain cross-contact cancellation before native storage. Unaffected
bodies preserve their exact native state and have no impulse receipt; a graph
involving only finite bodies leaves the rod unchanged. An all-zero solution
performs no native impulse event.

Publication checks every actual directional endpoint velocity, the actual
stored linear and angular transfer of each affected body, aggregate world
momentum including original support reactions, and cancellation-stable kinetic
work against both the contact work and its admitted loss. Opposite unresolved
body transfers cannot cancel into an accepted global account. All participants
return together only after those checks and a final cancellation checkpoint.
Poses, constitutive history, clocks and integration step counts are retained;
the returned native states continue through their original integrators. The
existing rod support restrictions remain binding.

`RodRigidGraphBudget` bounds body input/output storage before graph allocation:
the default and hard maximum are 128 bodies, including unused participants;
zero admits rod-only self-contact. The original normal law bounds the graph to
64 scalar contacts. Assembly needs O(contacts squared) storage/work and
O(contacts) native inverse-tensor actions; native staging needs
O(rod nodes + bodies + contacts) work without per-edge trajectory clones.
Collision discovery, event timing, persistent friction and distribution of
mechanical loss to heat or damage remain with their existing owners.

G0/G1/G3/G4 tests compare the singleton with the original isolated owner, a
mixed four-participant chain with independent simultaneous collision equations,
and eccentric multi-body attachments with an independent six-dimensional
material-origin spatial-mass solve. They check full tensor response, world-frame
and body/contact permutation covariance, native continuation, distinct-node
self-contact, rigid-only redundant cycles, support force/couple closure,
per-body unresolved-transfer refusal, bounded input, invalid geometry,
cancelled/failed publication and deterministic retry. A sequential-subgroup
falsifier establishes why the complete signed contact graph is required.

## Simultaneous Coulomb contact graphs on native rods and finite bodies

`rod_motion::RodMotion::apply_rigid_frictional_contact_graph` resolves normal
and two-tangent impulses across an entire graph containing the original rod and
multiple finite rigid bodies. Edges can connect a rod attachment and a rigid
body, two distinct rigid bodies, or two distinct rod nodes. The original
positive/negative endpoint identities determine the signs of both velocity
observations and reciprocal forces. Contact between two attachments on the same
rigid participant remains inadmissible.

`rigid_impact::graph::friction::RodRigidFrictionalGraphContact` adds a caller
tangent direction and an event-local resolved Coulomb coefficient to the common
graph geometry. The published `RodRigidFrictionalControls` and
`RodRigidFrictionalError` remain the controls and error boundary. Shared graph
preparation supplies the original endpoint geometry and compensated relative
velocities. The original frame owner supplies each work-conjugate contact
frame. Shared signed graph mobility retains full original tensor inertia,
eccentric arms, supported coordinates, all shared-participant cross terms, and
independent tangent frames at different edges.

The complete contact-major, normal-first `3N × 3N` matrix enters the original
`fs_contact::frictional_multicontact::solve_frictional_contacts` owner. Its
non-associated Coulomb law, Newton restitution, bounded nonlinear solve, global
passivity, determinacy and total-work admission remain authoritative. Initially
separating geometrically active contacts participate in that same solve.
Redundant impulse coordinates retain the original `UniqueVelocity` and
range/work witnesses; they are a representative, not unique local contact loads,
pressures or per-edge heat allocations. Original unresolved/nonunique outcomes
remain refusals.

Common graph staging applies one original rod point-impulse batch and one
aggregate native momentum update per affected rigid body. It checks the actual
post-event velocity triple at every edge, actual stored momentum increments,
stationary support reactions, aggregate kinetic work and the contact owner's
original relative work and loss before returning any participant. Tangential
impulses require actual witness collocation within the declared geometry bound.
Unused bodies retain their state and input order. Rigid-only graphs leave the
rod untouched; rod-only self-contact admits an empty body list. Zero-impulse
graphs invoke no native impulse owner.

The wrapper checks the original 21-edge / 63-row law bounds and the common
explicit participant budget before allocating graph arrays. Shared assembly is
quadratic in edge count with linear original inverse-inertia actions. Native
staging is linear in rod nodes, body count and edges, without a trajectory clone
per edge. Cancellation and any failed geometry, constitutive or balance check
return no candidate. Poses, material history, physical clocks and accepted-step
counts remain those of the original participants, which continue on their
native integrators. Collision detection, event timing, moving supports,
persistent friction history, multiple independent rod owners and thermal
partition remain separate contracts.

Focused G0/G1/G3/G4 groups compare a mixed three-mass chain to analytic
simultaneous stick/slip and falsify sequential pair updates; reconstruct
eccentric three-participant response from independent six-dimensional origin
mass blocks in sticking and sliding regimes; exercise world/tangent-frame,
edge and body-order changes; verify rod-only and rigid-only limits, original
support reactions, redundant-load determinacy and original refusal parity; and
check native continuation, bounded work, cancellation, exact no-ops and atomic
late refusal/retry.

## Simultaneous contact networks across independent native rods

`rod_motion::rigid_impact::graph::network::apply_rigid_normal_contact_network`
resolves one original normal-contact event across a bounded slice of continuing
`RodMotion<I>` trajectories and original finite `RodRigidParticipant` bodies.
`RodRigidNetworkEndpoint::Rod` identifies both a rod and its material node;
matching node indices on different rods remain distinct kinetic blocks.
`RodRigidNetworkContact` admits rod–rod, distinct-node self-contact, rod–body,
and body–body edges in the same original simultaneous law. Two attachments on
the same rigid participant remain inadmissible.

The rods share their inertia representation type, but keep independent native
models, section/material sources, masses, full central tensors, centroid
offsets, stationary supports, poses and integration histories. Rod attachment
offsets remain relative to the original material node, while rigid-body arms
remain COM-relative. Every rod referenced by a witness must represent the same
physical event time. Unreferenced rods preserve their own clocks without
constraining unrelated contacts. No pose, clock, constitutive history or
accepted-step count advances at this frozen-configuration event.

One shared graph core owns the signed original `J M^-1 J^T`, original point
kinematics, native impulse staging and actual-state admission for both indexed
networks and the existing one-rod graph APIs. Shared kinetic blocks contribute
their correct endpoint signs and full translation/rotation response; distinct
rod owners contribute no artificial cross mass. The original normal solver
retains complementarity, globally admissible restitution, total work/passivity,
and unique-impulse versus unique-velocity determinacy. Redundant impulse values
remain representatives, not unique local pressures or heat allocations.

Each affected rod receives one original point-impulse batch, and each affected
finite body receives one original resultant update. Every measured directional
contact velocity must match the admitted law. A single double-double account
then combines actual stored increments from all rods and bodies, including
each original rod's support force/couple and each body's independent requested
transfer. Rounded per-rod energies are never added to certify collective work.
The original native impulse gates remain authoritative. Common momentum cannot
inflate a transfer or work allowance, and failed numerical resolution cannot
be reported as physical dissipation. All native candidates return together
only after the law, endpoint, transfer, support-inclusive momentum, collective
kinetic-work and final cancellation checks succeed.

`RodRigidNetworkTransfer` preserves per-rod and per-body receipt ordering;
untouched participants have no impulse receipt. Zero solutions invoke no
native impulse owner. Empty rod/body slices are admitted when no witness names
them. The existing single-rod graph APIs delegate to this same owner and retain
their report shapes, original capacity and native behavior.

`RodRigidNetworkBudget` checks input/output participants before graph allocation:
the default and hard ceilings are 128 rods and 128 finite bodies, separately.
The additional configurable total-node bound includes untouched rods because
their continuing motion is returned; its default is 1,048,576 nodes and checked
addition prevents overflow. The original normal law remains bounded to 64 rows.
Matrix assembly uses quadratic contact work/storage and a linear number of
original inverse-inertia actions. Native staging is linear in total rod nodes,
bodies and contacts, with no per-edge trajectory clones or per-rod contact
solver. Collision discovery, event timing, moving supports, continuous contact
force evolution, heterogeneous inertia representation types in one slice, and
thermal/material-loss partition remain separate contracts.

G0/G1/G3/G4 groups compare three independently modeled rods and a finite body
with analytic simultaneous collision equations for elastic, partially elastic
and inelastic events; falsify sequential subgroup impacts; and check rod/body/
contact permutations and native continuation. Independent six-dimensional
material-origin spatial-mass blocks test eccentric rod–rod/body exchange with
full tensors under rotated, translated and boosted frames. Supported rods
independently close their actual force and couple reactions against one shared
body. Additional cases retain a known small collision loss during large common
momentum, preserve original single-rod and true rigid-only behavior, retain
redundant-cycle velocity determinacy, and exercise clock/participant/node
admission, cancellation, exact no-ops and deterministic retry.

## Simultaneous Coulomb contact networks across independent native rods

`rod_motion::rigid_impact::graph::network::friction::apply_rigid_frictional_contact_network`
resolves one simultaneous normal-and-two-tangent event across a slice of
independent native rods and a slice of original finite rigid bodies. A
`RodRigidFrictionalNetworkContact` adds a tangent hint and an explicitly resolved
event-local Coulomb coefficient to the indexed `RodRigidNetworkContact`.
Rod/node indices select an original trajectory and its original kinetic block;
equal node numbers or equal material/inertia values on different rods never
merge participants. Rods share their inertia representation type but retain
independent models, material properties, support sets and integration histories.

Cross-rod contact, distinct-node self-contact, rod/body and body/body edges may
coexist in the same original Signorini--Coulomb solve. The shared network owner
supplies actual material-point geometry, compensated relative velocities and
the complete signed contact-major `3N x 3N` mobility. Original full tensors,
eccentric offsets, stationary supports and shared-participant cross terms enter
without surrogate bodies or pairwise resolution. The published single-rod
`RodMotion::apply_rigid_frictional_contact_graph` now delegates to this owner;
its public controls, result and original single-rod capacity remain intact.

Every rod referenced by a witness must represent the same physical event time.
Unreferenced rods impose no clock constraint and retain their original state.
The existing contact owner retains its non-associated friction law, Newton
restitution, global passivity, determinacy, work witnesses and refusal classes.
Redundant impulse representatives do not determine unique local contact loads
or per-edge heat. Tangential impulses require actual witness collocation within
the declared geometry allowance; no endpoint is relocated to conceal a couple.

Common staging executes one native point-impulse batch per affected rod and one
native resultant momentum update per affected finite body. Before returning any
participant it checks the actual endpoint velocity triples, each original
transfer, complete stationary-support reactions, aggregate linear and angular
momentum and collective kinetic work against the original contact work/loss.
All continuing rods and bodies return in their input order. Poses, constitutive
state, clocks and accepted-step counts remain unchanged. A zero-impulse result
invokes no native impulse owner; cancellation or any early/late refusal returns
no partial candidate.

The original 21-contact / 63-directional-row solver limits, explicit rod/body
caps and checked aggregate rod-node bound are admitted before graph arrays are
allocated. Shared response assembly is quadratic in contacts with a linear
number of original inverse-inertia actions. Native staging is linear in total
rod nodes, bodies and edges, with no per-edge trajectory clone. Empty rod or
body slices are supported when no witness references them. Collision discovery,
event-time location, persistent friction, moving-support work and local thermal
partition remain separate owners.

G0/G1/G3/G4 groups compare a two-rod/body chain with independent simultaneous
sticking/sliding equations over restitution and friction regimes, and falsify
sequential pair resolution. An independent six-dimensional material-origin
mass oracle covers two eccentric rods, a full-tensor finite body, sticking and
sliding, world/tangent-frame changes and rod/contact permutations. Further
fixtures reconstruct cross-rod support force/couple reactions; retain original
redundant-load determinacy and constitutive ambiguity refusals; and check native
continuation, participant identity, clock compatibility, bounded nodes/owners,
late imbalance refusal, cancellation/retry, exact no-ops, a true rigid-only
empty-rod call and published one-rod API equivalence.

## Coulomb friction in isolated native rod/rigid impacts

`rod_motion::RodMotion::apply_rigid_frictional_impact` resolves the normal and
both tangent impulses together through the existing full three-dimensional
contact law. Its published `RodRigidFrictionalContact`,
`RodRigidFrictionalControls`, `RodRigidFrictionalReport` and
`RodRigidFrictionalError` API remains intact. The explicit coefficient is
resolved by the caller for this instantaneous event. Original contact-frame
admission supplies the same axes for observations and reciprocal impulses.

The complete symmetric three-by-three mobility retains eccentric arms, full
tensors, normal/tangent cross terms and supported rod coordinates. The report
retains separate original rod and rigid mobilities and the actual world
impulse. Nonzero tangential impulses require collocated original witnesses
within the declared geometry tolerance; no midpoint relocates either source.
The actual relative work is checked independently against the original contact
owner's work as well as actual native kinetic change and total impact loss.

Normal, frictional and aggregate rod/body events now share original point
preparation and an all-node native momentum/work audit. Compensated relative
translation is formed once before projecting onto contact axes. The audit
checks actual stored linear and angular increments, stationary support forces
and couples, measured endpoint velocities and polarized kinetic work. It keeps
small loss against large common drift and refuses an unresolved native impulse
even when an angular response remains representable.

Both continuing states publish after final cancellation and balance checks.
No-op states retain their bits; accepted events retain original geometry,
material history, clocks, step count and compensated drift remainder.
Constitutive iteration, passivity and determinacy refusals remain with
`fs_contact::frictional_multicontact::solve_frictional_contacts`. Total
mechanical loss does not imply a thermal partition, persistent friction history,
coefficient evolution, collision discovery or contact duration.

The original four friction test groups remain intact. Additional independent
spatial-mass sticking, large-common-drift loss and rotational-approach refusal
fixtures exercise the shared preparation and audit together with the
simultaneous normal-contact extension.

## Simultaneous native rod and finite rigid-body normal impacts

`rod_motion::RodMotion::apply_rigid_normal_impacts` resolves a localized set of
normal contacts between one original native rod and one finite rigid body.
Each witness keeps its own original material-node attachment, body-frame rigid
arm, normal, proximity and restitution. The complete Delassus matrix includes
the shared finite body's translation and full tensor rotation. Contacts on the
same rod node also share its admitted translational coordinates and full
central tensor response, including eccentric material-node offsets. Different
rod nodes retain their original independent instantaneous inertia.

Inverse-inertia actions are prepared once per contact direction. Matrix
assembly then has quadratic contact-sized work and storage, without a full rod
matrix or repeated factorization. The adapter checks its 64-contact limit and
the original caller variable cap before allocation. Cancellation is polled
through preparation, assembly, native staging and publication.

The existing `fs_contact::multicontact::solve_normal_impacts` owns the globally
coupled normal complementarity and declared restitution. Initially separating
but geometrically active rows remain in the solve; another row can make them
close. Restitution choices that inject assembled kinetic energy refuse. Contact
ordering does not replace the simultaneous law with a sequential collision
scheme. The original PSD determinacy witness is retained: a `UniqueVelocity`
result has a representative impulse vector, not unique local contact loads or
heat allocations.

All selected impulses act through one original `apply_point_impulses` call.
Their equal and opposite world actions and body-frame moments accumulate with
compensation into one native rigid `apply_momentum_impulse` call. The moment
uses the same rounded physical impulse passed to the rod. The shared native
audit checks actual stored changes across all affected rod nodes, original
stationary support reactions, finite-body linear/angular momentum, every
post-event relative normal speed, and compensated pair kinetic change against
relative impulse work and the original contact loss.

Both continuing states publish only after final cancellation and all physical
acceptance bounds pass. Poses, material history, physical clocks, accepted-step
counts and compensated drift state are preserved. An empty or all-zero event
still admits the original controls and returns unchanged states without native
impulse receipts. No collision discovery, event timing, penetration repair,
friction law, multiple rigid participants or thermal loss destination is
inferred.

`tests/rod_motion/rigid_impact/multiple.rs` checks analytic two-attachment
collisions at plastic, partial and elastic restitution; falsifies sequential
pair treatment; checks initially separating activation and contact permutation;
compares eccentric same-node response to an independent six-by-six spatial
mass solve with full rigid tensors and rotated/boosted geometry; retains
redundant-load determinacy and stationary support reactions; and exercises
native continuation, passivity refusal, budgets, cancellation, retry and no-ops.
This extends the native flexible-body transfer in `frankensim-b8bxd.9.3` while
keeping the existing contact owner responsible for the constitutive law.

## Simultaneous rod and finite rigid-body Coulomb impacts

`rod_motion::rigid_impact::friction_multiple` extends the original isolated
frictional impact to a set sharing the same native `RodMotion` and the same
finite `RigidBodyState`. `apply_rigid_frictional_impacts` consumes the published
`RodRigidFrictionalContact`, `RodRigidFrictionalControls` and
`RodRigidFrictionalError` contracts. The complete contact-major, normal-first
`3N × 3N` Delassus matrix includes the original rod's full central tensor and
eccentric material-point arms, the supported translational/rotational freedoms,
same-node cross response, and every cross response through the shared rigid
body's full tensor. Distinct rod nodes share the finite body; their nodal mass
blocks do not acquire fictitious direct coupling during the instantaneous event.

`fs-contact::frictional_multicontact` retains ownership of the non-associated
Signorini–Coulomb law, Newton restitution, the bounded nonlinear solve, global
passivity, numerical determinacy and total-work admission. All geometrically
active witnesses participate, including initially separating contacts. The
adapter returns the original solver receipt and complete original matrix;
`UniqueVelocity` retains its range/work witnesses and does not authorize a
unique local force, pressure or per-contact heat allocation. Unsupported,
unresolved or nonunique original responses retain their original refusal.

Accepted world point impulses enter one original rod point-impulse batch.
Opposite body forces and central couples are compensated-summed and enter one
original rigid momentum-impulse update. Each actual normal and tangent endpoint
speed is reconstructed from the stored candidate states. The common native
audit checks support-inclusive world linear/angular momentum, actual kinetic
work, the original contact relative work and total impact dissipation before
either candidate is returned. Any nonzero tangent impulse requires the
original witness points to be collocated within the declared geometry bound;
their actual separation still enters momentum admission. Stationary support
reactions are retained and do zero work.

The original contact cap is 21 witnesses / 63 scalar rows, further limited by
the caller's contact and matrix budgets before allocation. Preparation uses
`O(N²)` work/storage and `O(N)` prelowered inverse-inertia actions. Publication
uses one native update per participant and `O(rod nodes + N)` rod work; it does
not perform per-contact state copies, elastic/source evaluation or a new mass
factorization. An all-zero solution preserves both native states without
calling their impulse owners. Cancellation and any failed law, geometry,
endpoint or balance check return no candidate. The event leaves configuration,
material history, time and step count unchanged, and accepted states continue
through the original time integrators. Collision search, gap repair, persistent
friction, thermal partition and contact sets involving multiple finite rigid
bodies remain separate contracts.

The original rod point-impulse owner now keeps compensated per-node force,
material-point moment and free-couple sums through the centroid shift and one
final native momentum rounding. All prescribed, normal-contact and frictional
batches share this behavior. Global momentum and work receipts use compensated
reductions as well. Actual stored increments remain subject to the original
per-node and aggregate transfer bounds; compensation does not admit unresolved
kicks. Focused original-owner fixtures retain a small net force/couple across
all six orderings of large cancelling impulses, verify eccentric mass and
support reactions, and retain small physical work after large opposing nodal
kinetic transfers. Scratch remains O(nodes), with O(nodes + impulses) work.

Acceptance coverage includes isolated-API equivalence, analytic coupled
stick/slip, an independent dense six-dimensional origin spatial-mass oracle,
same-node and distinct-node contacts, full eccentric tensors, permutation and
world/tangent-frame changes, supported and redundant contact sets, preservation
of the original unresolved/nonunique refusal, initially separating and inactive
witness response, native continuation, bounded work, cancellation, exact no-ops
and atomic late refusal/retry.

## Source-backed ellipsoidal fluid carriers

`material_immersed::RetainedFluidEllipsoids` binds one centered uniform
source-resolved region per carrier to `VolumeMoments::ellipsoid` geometry. An
explicit proper rotation maps principal shape axes into material/body
coordinates; sampled geometry and the complete source-derived inertia tensor
use that same rotation. The shared SI owner advances these carriers through the
existing simultaneous D3Q19 moving-body/full-tensor solver, retaining source
receipts, atomic refusal and same-density viscosity rebind. This remains
retained-interior-fluid coupling plus additional carrier inertia. Native
enclosing-sphere separation and wall/halo admission apply; spherical lubrication
and a displaced-fluid closure are not inferred. The original sphere API shares
the same owner and retains its original film path. Tests cover analytical rotated
inertia, frame-equivalent populations and motion, SI fluid/body angular balance
including periodic transport, and state preservation under rebind and refusal.

## Native circuit dielectric absorption and polarization heat

`material_circuit::dielectric::DebyeCircuit` replaces one original capacitor
placeholder with its source-derived C_inf and appends the source's series RC
polarization arms without changing original branch indices. Explicit original
capacitor charges/inductor flux linkages and separate polarization charges enter
the existing `CircuitDae` consistent initialization unchanged. The prepared
descriptor is shared across candidates; every step delegates to the original
audited solver rather than adding a dielectric time integrator.

The source-bound polarization charges carry memory, so a released prepolarized
capacitor recovers terminal voltage from conserved total charge. Native midpoint
voltage drops independently determine each process's positive heat. That heat
is already inside total circuit dissipation and is reported separately from
other authored resistor losses. Both the dielectric's terminal-work/storage/heat
account and the complete circuit balance must close within explicit absolute
budgets. Receiving thermal owners must accept the returned heat with the same
electrical candidate; it is not fabricated bulk conductivity or automatic
finite-temperature feedback.

Storage increments use compensated charge/flux changes times midpoint efforts,
including the midpoint rounding correction; endpoint energies remain diagnostics.
The dielectric checks its own storage/work/heat account, and an independent sum
across native capacitors and inductors checks the descriptor's cancellation-stable
storage increment. Small polarization heat and inductor work remain observable
below the spacing of a large pre-existing electric energy. The bounded
polarization endpoint buffer is allocated on the stack.

The complete original dielectric material bundle and fixed geometry remain
attached. Original other R/L/source/transformer branches retain their authored
values; the ordinary `MaterialCircuit` conductor-only resistor binding is not
relabelled to accept polarization friction. Nonphysical parameters, bounded
graph/state overflows, inconsistent ideal-source impulses and cancellation
refuse before publication. Time accuracy needs relaxation-scale resolution and
step refinement independently of algebraic energy checks. This base isothermal
API supplies no switching impulse, finite-temperature feedback, variable
geometry, saturation or breakdown; finite thermal feedback is provided by the
separate `circuit_enthalpy::dielectric` owner described below.
G1/G3/G4 `tests/dielectric_circuit.rs` checks voltage-clamped multipole relaxation,
open-circuit absorption with charge conservation, second-order refinement,
current-driven frequency/time agreement and loss, the zero-memory native RLC
limit, separate load heat, small polarization and inductor work under a large
held bias, and atomic refusal. All materials are synthetic.

## Simultaneous rigid and finite-solid frictional impacts

`finite_solid::contact` connects the original solid impulse/mobility interface
to native rigid-body states and `fs-contact::frictional_multicontact`. One solid
can interact with multiple rigid bodies in a single globally admitted event.
The full contact matrix retains solid cross-point coupling and rigid full-tensor
translation/rotation response, including normal/tangent cross terms. All normal
and Coulomb impulses are solved together by the existing contact owner, with its
original determinacy and passivity refusals. A velocity-only result retains
representative, non-authoritative individual impulses and unpartitioned loss.

`FiniteSolidBody::rigid_contact_candidate` applies equal/opposite impulses at the
original body arms and current barycentric solid points. The solid receives the
complete batch in one native mass solve. Actual endpoint velocities, stable
kinetic increments, world spin/orbital momentum, stationary-support reactions
and joint work are checked before either candidate is returned. An orthonormal
normal-first frame is explicitly admitted. Inputs remain unchanged on source,
geometry, solver, physical-balance or cancellation failure; poses and time do
not advance. Callers provide event geometry and interface coefficients. No
collision detection, penetration repair, continuous contact trajectory or loss
to heat conversion is inferred. A finite activated gap is retained in the angular
balance and can reject a tangential couple; witnesses are never relocated.

Controls bound rigid bodies, contacts (through the original 21-contact owner),
each consistent-mass solve and independent physical acceptance tolerances.
Tests exercise source-mass stick/slide/restitution, two rigid bodies coupled
through negative solid cross mobility, off-center full-inertia impact, fixed
support reactions, canonical contact ordering and late rejection of a real
gap-induced angular couple. They use independent tetrahedral mass and impulse
oracles; native execution remains unverified in the current resource-blocked
environment.


## Harmonic Rayleigh pressure and vector receivers

`pcm_wav::baffled::harmonic::Prepared` is the continuous-frequency image of the
existing `BaffledPressure` surface basis. It consumes the same positive physical
areas and complete signed `SurfaceSample::mode_shape` coefficients in the z=0
infinite rigid baffle. A 1/area basis still denotes volume flow. Prepared transfer
columns expose pressure, Cartesian pressure gradient and fluid particle velocity
for every input and receiver; repeated coherent source evaluations allocate no
memory. The shared `fs_bem::near_field::FirstOrder` produces ideal omni, cardioid
and figure-eight Pa-equivalent observations with the same resolved medium.

All phasors use `exp(-i omega t)`. For one volume-velocity element Q,
`p = -i rho omega Q exp(ikR)/(2 pi R)` and
`u = Q exp(ikR) (1/R^2 - i k/R) r_hat/(2 pi)`.
The gradient is `i rho omega u`. Retaining the velocity's reactive `1/R^2` term
gives `u_r/(p/(rho c)) = 1 + i/(kR)`, so proximity response and phase follow from
the field itself. Signed complex surface/input contributions are summed before
microphone projection; a scalar pressure null may retain transverse particle
velocity and directional output.

Preparation admits at most the original 128 input coordinates, the shared BEM
64 receivers and 120,000 surface points. Complete sample/input/receiver work is
checked against the declared budget (default 2^21, hard cap 2^24 contributions),
with cancellation before receiver rows and every 64 surface samples. The default
receiver clearance is the existing 5 cm; an explicitly smaller positive clearance
requires finer caller-supplied quadrature. Invalid/nonfinite inputs, unrepresentable
fields, clearance violations and incomplete arrays refuse. Source evaluation and
microphone projection leave output buffers unchanged on refusal.

Focused native tests compare the analytic baffled monopole and ideal cardioid
for kR=0.05..50, the far-field velocity limit, complex spatial finite differences,
a pressure-cancelled but velocity-live source pair, the original circular surface
quadrature against the analytic piston, and the original differentiating causal
pressure owner through 21.6 kHz at 48 kHz (relative RMS difference <5e-4 after its
reported numerical latency). Full 64-receiver/128-input order, work and cancellation
checks use the same production kernels and original `FirstOrder` implementation.

This provides exact Green kernels at the supplied quadrature points, not an error
bound for the continuous surface. Structural and spatial refinement remain the
caller's responsibility. The harmonic matrix has no sampling clock or temporal
interpolation filter. It adds no causal vector receiver, measured capsule response,
finite-baffle/room scattering or reciprocal radiation load, and does not close
the wider M2.6-W1 work item.

`Prepared::new_relaxing` additionally binds the original evaluated
`fs_material::relaxation::RelaxationState`, deriving density and reference sound
speed from that same source. Its `causal_wavenumber` supplies molecular and
classical transport loss/dispersion once per frequency, shared with the existing
time-domain gas observer. Each original surface-to-receiver path uses its own
`exp(i k R)/R` before coherent summation; the velocity factor is the complete
`1/R - i k`, including absorption's spatial derivative. The immutable source and
used complex wavenumber remain inspectable. Existing work/input/receiver limits,
cancellation, transactional buffers and allocation-free warm evaluation apply.
The explicit lossless constructor retains its original arithmetic.

This is the original effective homogeneous gas-propagation image, with real
equilibrium density in its Rayleigh source and momentum relation. It does not
claim full thermoviscous boundary-layer fields, universal high-frequency
transport validity or equivalence to `fs-bem`'s reservoir-only molecular reaction
model. Ideal microphone velocity sensitivity retains the source's equilibrium
`rho*c` reference. G1/G3 integration tests use an independently evaluated synthetic
relaxation polynomial, outward intensity decay, signed unequal propagation paths,
spatial finite differences, and the original causal material-pressure observer.

`Prepared::{new_thermoviscous,new_thermoviscous_aperture}` additionally accept
the original `ClassicalBulkAbsorption` state of a liquid or real gas. Its
`thermoviscous_wave` owns both the complex wavenumber and thermodynamic-pressure
impedance `Z_p`; no second absorption model or ideal-gas EOS is introduced.
For the prescribed-flow potential `phi=-Q exp(i k R)/(2 pi R)`, the shared
quadrature now retains `u=grad(phi)`, `p=i k Z_p phi` and
`grad(p)=i k Z_p u`. Replacing `k Z_p` by `rho omega` would confuse pressure
with longitudinal traction in a viscous medium. Lossless and effective molecular
paths retain their prior arithmetic. The common `FirstOrder` receiver still
weights velocity by real equilibrium `rho*c`; strong loss does not retain the
lossless plane-wave cardioid normalization. Point and finite-area evaluation
share the same budgets, cancellation and atomic output publication.

Five new focused native tests cover independent Cp=Cv pressure/velocity limits
from dimensionless viscous loss 1e-8 through 1e12, original viscous momentum and
spatial derivative refinement, a four-rung disk ladder against closed complex
area integrals, lossless convergence, singleton identity and transactional
refusal. Six original harmonic regression tests also pass through the same
source numerics. This focused native closure is not a full Cargo workspace gate.
The result is the homogeneous longitudinal acoustic branch on the ideal baffle;
independent thermal/shear modes, no-slip/thermal boundary layers, material
frequency support, finite-baffle scattering and backreaction remain separate.
It does not add a causal material-directional receiver or close M2.6-W1.

### Finite-area coherent microphone capsules

`pcm_wav::baffled::aperture::Aperture` describes positive sensing-area
quadrature in the original source frame. Its circular builder exposes radial
and azimuthal refinement independently. The measured quantity is the area
average of pressure and particle velocity, followed by the original common
`FirstOrder` pattern; it is never an average of magnitudes, powers, source
distances or independent recordings. The plane normal and microphone front
axis are separately declared. The actual source-to-capsule paths retain their
own phase, travel time and reactive velocity, so a pressure null may retain
directional output. DD-normalized area weights survive until the complete
coherent field or FIR coefficient is rounded.

`BaffledPressure::{new_aperture_first_order,new_multirate_aperture_first_order}`
compile the complete surface/capsule product into the original one pressure
bank and optional aligned reactive-velocity bank. Their callback and physical
clock are unchanged. The reported position is the area centroid for metadata;
it never replaces the physical path geometry. Every point retains the 5 cm
baffle-clearance admission. Up to 4096 sensing points are accepted, with an
explicit complete coefficient-projection work limit, bounded history, and
cancellation checks before/after point-template preparation and during merged
path projection. A singleton delegates bit for bit to the point constructor.

`harmonic::Prepared::{new_aperture,new_relaxing_aperture}` tile the existing
Green-kernel owner in groups of at most 64 receivers, then retain just one
pressure/vector row per source input. Warm source evaluation allocates nothing
and is independent of capsule quadrature count. Molecular propagation uses the
original source-owned complex wavenumber on each path, including its gradient;
no second absorption law is applied. The harmonic work budget covers all
source/input/capsule contributions before any tile is prepared.

Focused G1/G2/G3/G4 tests cover a three-rung disk ladder against independent
closed-form pressure and axial-velocity area integrals, coherent pressure
cancellation with live velocity, a 2^-60 sensing-area residual, singleton
equivalence, multirate waveforms against the continuous vector field, and
complete-work/clearance/cancellation refusal. Source and capsule spatial
refinement remain separate caller obligations. These are ideal stationary
sensing apertures in the same infinite rigid baffle model, not diaphragm
mechanics, capsule-body scattering, measured capsule frequency responses,
acoustic backreaction or closure of the broader M2.6-W1 bead.

## Three-dimensional principal-axis moisture transport and swelling

`moisture::cartesian::MoistureCartesian` evolves dry-basis bound moisture on
actual uniform Cartesian reference volumes. Its three explicitly assigned
source diffusivities are the diagonal principal values aligned with the named
grid frame. Repeating a source key explicitly assigns the same value to multiple
axes; no isotropy or rotated tensor is inferred. The shared slab/grid binding
admits typed absolute temperature and dry-mass moisture, pins original material
claims, and obtains dry density only from the reference wet density divided by
`1 + reference_moisture`. Dry skeleton mass stays fixed.

By default, internal faces use symmetric harmonic series conductance. External
faces, in -x,+x,-y,+y,-z,+z order, use the original prescribed-moisture, sorption
and finite surface-emission laws, including half-cell diffusion resistance.
Each cell stores water once even at edges/corners or when an axis has one cell.
The seven-point matrix and Jacobi preconditioner require O(N) storage and work
per application; the original checked-CG owner performs each bounded Picard
solve. Endpoint source refresh and an independent Euclidean equation residual
must pass alongside the water budget computed from those same endpoint face
fluxes. All source/domain, maximum-principle, increment, work and cancellation
refusals leave the accepted field unchanged. Actual interface moisture is also
re-admitted by the original source. Reports retain cell and surface source
states, six physical transfers and native linear/nonlinear work separately.

`MoistureSolid::from_cartesian`, `resolve_cartesian` and `advance_cartesian`
map each original moisture cell to six conforming tetrahedra in the complete
3D lattice. They use the existing local elasticity, integrated directional
swelling, compatible displacement, stress, geometry and prestressed wet-inertia
owners without averaging through thickness or across width. The grid and solid
mesh frame must agree. A mechanical refusal discards the combined endpoint.

Focused authored regressions cover a separable three-axis anisotropic Neumann
eigenmode, all six single-cell boundary transfers, finite-emission and
impermeable source faces, nonlinear endpoint diffusion
and bounded refusal, reduction to the original slab, and through-thickness
uptake driving local native swelling and compatible shape. These are numerical
and source-contract tests, not measured-wood validation. The time method remains
first-order backward Euler; CG convergence is not a spatial/time error estimate.
The spatial owner remains homogeneous, isothermal and fixed-reference, with
principal axes aligned to the grid. It does not add unlike-material partition
laws, arbitrary tensor rotation, sorption heat,
stress-assisted diffusion, hysteresis or finite-deformation transport, and does
not close the broader MR46 measured-material acceptance item.

### Nonlinear principal-axis moisture potentials

`MoistureCartesian::new_kirchhoff` and `with_kirchhoff_diffusivity` extend the
homogeneous source-potential scheme to the three explicitly authored,
grid-aligned principal laws D_i(M). Each trial prepares one original material
chart per axis over its complete maximum-principle moisture interval. Internal
fluxes use exact directional potential secants; each physical boundary reuses
the bounded nonlinear half-cell/surface-emission balance. Charts are reused
across Picard iterations without per-face source resolution or query-map
allocation, and the report carries their x/y/z identities.

Nonnegative laws may wet a dry D_i(0)=0 state; identically zero directions remain
impermeable. The original symmetric seven-point operator, native checked CG,
nonlinear/residual/water/maximum-principle gates, cancellation checkpoints and
transactional publication remain authoritative. Existing finite-reservoir and
solid consumers retain the selected field mode. This does not supply rotated
or spatially heterogeneous tensors, unlike-material partitioning, or a
continuum/front-location certificate.

Five authored physical tests cover independent six-face anisotropic uptake,
simultaneous three-direction interior dry-corner wetting, zero-transverse
reduction to native slabs, six nonlinear finite-emission interfaces, and
bounded refusal without state publication.

Finite reservoirs consume this same Cartesian owner through
`FiniteMoistureAtmosphere::advance_cartesian` / `advance_shared_cartesian` and
`FiniteMoistureSolution::advance_cartesian`. A private field adapter resolves
the original face sorption/finite-emission laws and invokes native slab or 3D
transport; it does not reimplement either reservoir's scalar controller. Public
candidate/report defaults preserve existing slab callers. Every reservoir trial
restarts the original spatial field, with six faces in -x,+x,-y,+y,-z,+z order
and complete endpoint nonlinear-source checks. Shared 3D specimens retain stable
ID ordering, common accepted clocks/temperature and an aggregate cell budget.
Local water increments debit the one vapor inventory once; original gas
reconstruction independently checks total mass and each noncondensing carrier.
Rejected bracket candidates release their full source fields; only residual
signs survive into bisection.

The solution controller retains fixed nonvolatile solute amounts while water
uptake concentrates and desorption dilutes the original solution. Its existing
pure-solvent path and bounded safeguarded-volume solve remain unchanged. This
is the equal-pressure, negligible-vapor-storage activity bridge, not liquid
immersion permeability. The Cartesian solid's `advance_cartesian_finite_atmosphere`
and `advance_cartesian_finite_solution` admit reservoir, spatial diffusion and
native swelling/stress/wet inertia together; any late source or mechanical
failure publishes none of them. Gas internal-energy change excludes bound-water
caloric storage and sorption heat; no full thermostat or chemical-work balance
is added. Focused authored regressions cover finite-vapor two-layer gradients
and late refusal, six-face nonlinear uptake, simultaneous shared adsorption and
desorption with finite/zero surface exchange, and manufactured solution uptake
and desorption reaching native 3D swelling. These remain numerical/source
verification fixtures, not experimental moisture-material validation.

### Whole-duration adaptive moisture evolution

`moisture::adaptive` adds `MoistureSlab::advance_adaptive_to` and
`MoistureCartesian::advance_adaptive_to`, plus corresponding
`FiniteMoistureAtmosphere` / `FiniteMoistureSolution` slab and Cartesian
methods. One private controller compares a native full interval with two
native half intervals and accepts only the two-half endpoint. It independently
checks every local dry-basis content and, where present, RH/activity and actual
external water mass. The relative allowance scales changes from the accepted
state; large static water inventories cannot hide a small unresolved transfer.
Unequal representable half intervals retain their first-order error factor.
The estimate is local and asymptotic, not a global accuracy certificate.

Only accepted fine fluxes enter compensated physical-face accounts. A single
fixed whole-duration tolerance checks original-to-final local bound-water
change against those transfers and, for finite reservoirs, against external
water depletion. Finite vapor independently rechecks whole-window total mass
and each noncondensing carrier with the original gas inventory owner. These
budgets are never multiplied by the accepted step count. The original native
step budgets, source admission and source-selected endpoint diffusivity remain
active at every probe. The temporal work ceiling counts every full/half native
field or coupled-reservoir call, including rejected/failed probes; the original
bounded scalar, nonlinear and linear controls also apply within each call.
Only typed convergence or content-increment refusals and temporal error permit
shorter retries. Source/domain, conservation and cancellation errors remain
terminal. No failure returns even a successful internal physical prefix.

Reports retain the complete accepted transfer totals, original/final water and
RH/activity, work counts, step sizes and the last native source receipt, without
retaining an unbounded source-state history. Boundaries and source temperature
stay fixed during a call; environmental changes require explicit interval
breaks. Native mechanics consume the final spatial field as before. This adds
neither thermal/sorption energy coupling nor time-varying source extrapolation.
Authored regressions exercise continuous-time discrete diffusion decay and
refinement, Cartesian/slab reduction, finite-vapor two-store relaxation,
constant-activity finite-solution water loss, typed increment refinement,
domain/cancellation refusal and a late work failure with no partial result.

`FiniteMoistureAtmosphere` and `FiniteMoistureSolution` additionally expose
`advance_shared_adaptive_to` and `advance_shared_cartesian_adaptive_to`.
They use that same temporal loop around the original simultaneous reservoir
roots. Private snapshots retain separate stable-ID-tagged native field views,
source bindings, cell counts and dry masses; no cells are concatenated into a
surrogate geometry. Every specimen's componentwise moisture discrepancy and
the one reservoir's RH/activity and water mass must satisfy their original
temporal scales. Neither specimen averaging nor a zero net bath change can
hide opposing local uptake/desorption errors.

Accepted fine transfers accumulate separately for every original face of
every specimen. `SharedMoistureAdvanceReport::specimens` retains each ID's
whole-window water change, face transfers, residual and maximum accepted
moisture error. The common summary independently closes original-to-final
stored water against all accepted transfers and reservoir depletion, with
the same fixed window tolerance rather than a participant/step multiplier.
Finite vapor also rechecks every original carrier and total mass. Aggregate
face slots are scalar sums of corresponding local slots, not a common spatial
surface. The original fields and their last native transport/source receipts
are returned as `SharedMoistureSolid`; source history stays bounded.

All specimens must share the original clock and reservoir temperature; every
half/full endpoint rechecks the actual common clock, original binding and
stable IDs. Aggregate cells are bounded before cloning and on every original
root call. The shared temporal evaluation ceiling counts complete coupled
probes, including rejected attempts, and the existing scalar/nonlinear/linear
bounds remain active inside each probe. A late source/domain, inventory, work
or cancellation refusal discards the complete requested window. The caller's
fixed physical source/volume/RH brackets must cover its entire duration.
Existing single-field APIs and report types retain their behavior through the
same controller. Focused authored regressions cover opposing-field refinement
with stationary pure solvent, finite-vapor mean/difference relaxation through
both geometries, stable-ID replay, signed per-specimen accounts, and refusal
after a successfully computed private prefix. Native moisture execution remains
a separate validation gate; these fixtures add no caloric or experimental claim.

## Complete modal reflection admission

MatrixFirLine's wall-normalized pressure waves all use the same positive real
zc0. Net mean input power is (||a||² - ||b||²)/(2 zc0), so admission uses the
largest singular value of the complete reflection matrix, including coherent
higher-mode inputs. Every raw frequency sample and the complete realized FIR
on the 4N grid (including DC and Nyquist) must be contractive within the existing
passivity budget. The shared complex SVD must converge; no entry or operator is
clipped into admission. This finite-grid check does not certify unsampled
frequencies. The cancellable constructor polls source bins, transform/tap tiles,
and SVD boundaries, with the existing modal and 65536-tap ceilings.

The matrix ladder checks all-entry tail energy and all-entry common causal-prefix
change, as well as plane resonance movement. One existing Nyquist taper is
computed at the ladder ceiling and subsampled at every rung, keeping the
numerical transfer function fixed and avoiding repeated taper FFTs for each
matrix entry. Other scalar realization callers retain their original rung-local
taper arithmetic. A small plane response or tail does not qualify a changing or
active higher-mode FIR. The current precursor treatment can amplify an otherwise
passive near-lossless evanescent port; such full-matrix realizations are explicitly
refused. Physical matched-port witnesses cover the affected lossless resistance
fixtures and Bessel brass/long-cylinder entrances. Successful one-mode matrix
runtime controls and a proven passive coupled-FIR push test remain; these do not
claim qualification of the refused higher-mode physical realizations.

The new Rust regressions have not yet been executed in a native Rust toolchain;
independent numerical reproductions do not substitute for that build/test gate.


## Nonideal finite-solution activity in moisture transport

`moisture::solution::FiniteMoistureSolution` now retains the same shared
`AqueousSolution` as finite-volume membrane flow. Existing ideal inputs convert
directly. A declared virial model supplies its own endpoint water activity and
osmotic-pressure diagnostic; the existing native source-backed sorption law
receives `RH = a_w` on the equal-pressure, isothermal, negligible-vapor-storage
bridge. Retained solute moles and the shared thermodynamic model survive every
candidate. Concentration-support refusal occurs before either accepted water
inventory changes. Native diffusion, surface-exchange sources, solvent-plus-solid
water balance, cancellation and downstream moisture-dependent plate/solid
response retain their existing owners.

This extension does not add liquid-immersion permeability, volatile-solute
exchange, ionic transport into the solid, vapor-space storage, sorption heat or
stress feedback into solvent chemical potential. New native regression sources
manufacture a nonideal one-cell endpoint and verify the original sorption source
queries, then exercise unsupported-concentration refusal and a supported retry.

### Source-resolved directional rigid/solid impacts

`finite_solid::contact::elliptic` extends the same native rigid/solid event to
directional Coulomb friction. `rigid_elliptic_contact_candidate` delegates the
complete contact matrix to the existing `fs-contact` ellipse solver; its physical
impulses enter the original rigid-body and consistent solid-mass owners. Shared
assembly and acceptance retain cross-point and normal/tangent coupling, full
rigid inertia, stationary-support reaction, actual endpoint velocities, stable
kinetic increments and world spin/orbital momentum. The circular API retains
its default types and behavior. `with_friction` binds a new law without changing
either original material witness.

`rigid_source_elliptic_contact_candidate` reuses the existing ordered-interface
binding and `SourceThermalEllipticCoulomb::sample_preimpact_coefficients`.
Geometrically active contacts are queried once in canonical ID order, including
initially separating contacts that shared mobility can close. Each query carries
the actual native point velocity and exact normalized contact frame. The source
owner substitutes actual incoming slip speed, retains accepted accumulated path
and texture, and checks original claims, authority, domains and complete query
context. The caller supplies admitted temperature and optional pressure and
explicitly selects the static or kinetic coefficient role for the fixed event.
Inactive contacts query no source and retain their accepted history.

The returned source samples contain candidate states indexed by original input
position; callers adopt these together with both mechanical candidates. Source,
geometry, solver, cancellation or late balance failure changes no borrowed state.
Changing a contact frame requires the existing objective history transport.
Transformed numerical determinacy witnesses remain explicitly labeled; unique
velocity does not establish unique individual impulses or contact heat shares.
This is a preimpact-sampled instantaneous law. It introduces no pressure from
impulse, event duration, impact heat partition, wear, penetration repair or
classical adjoint.

Focused G1/G3/G4 tests cover analytical off-center source-speed coupling,
temperature/role/texture substitution, a rotated ellipse with full three-axis
mobility, active/inactive source ordering, original claim retention and source,
finite-gap angular-balance and cancellation refusal. Coefficients are synthetic.

## Permanent creep on the spatial solid/enthalpy history

`solid_enthalpy::creep::CreepSolid` binds one original Norton
material source and actual caloric cell to each native reference tetrahedron.
With an unwrapped Norton source, this owner selects zero thermal expansion,
constant elastic storage and fixed reference density/volume. The material owner's shared caloric
admission verifies the original card, complete fixed source context, pinned
Cp/density claims and wholly solid supported interval. No missing coefficient
is inferred from a material name or a melting point.

Each interval reuses `SmallStrainTetProblem::solve_creep_from` and the original
`DepositionSimulation`. Local temperature controls Arrhenius kinetics; native
stress redistribution and retained permanent creep change the actual loaded
shape. Every coupling trial starts from the same accepted displacement, tensor
history, enthalpies, observer and clock. Final mechanics are solved again at
the returned spatial temperature field. Per-cell absolute heat defects prevent
cancellation between underheated and overheated cells; separate temperature,
thermal and combined work gates precede publication of the complete candidate.

Only physical creep dissipation enters the existing thermal source hook.
Backward-Euler loss, applied/support work and elastic storage retain their
native distinct accounts. All existing observer hooks, including interval
sources and reservoir conductances, are forwarded; conduction, boundary heat,
events and additional thermal cells retain their original owner and energy
account. This is a quasistatic small-strain reduction; time/spatial refinement,
cell temperature uniformity and fixed thermal geometry remain caller/model
obligations. It adds no phase change, inertia, damage or service-life claim.

Manufactured `tests/solid_creep_enthalpy.rs` compares heated constant-stress
creep with an independent scalar thermal/kinetic equation, varies material rate
and thermal substeps, conducts heat through a nonuniform six-tetrahedron loaded
bar, checks original source identities, and exercises reservoir/interval-source
composition and transactional failure/cancellation.

## Retained nonlinear material-circuit capacitance

`MaterialCircuit::new_with_capacitor_curves` binds explicitly authored monotone
V(q) curves through the original native descriptor alongside source-resolved
bulk/contact resistance and signed mutual inductance. Original capacitor branch
indices, charges, flux linkages, ideal-source ports and transformer topology are
retained. Ordinary material, contact, conductivity and finite-enthalpy trials
keep the same prepared charge laws; they do not silently return to the graph's
linear placeholder capacitance. Immutable prepared knot/slopes/energy tables
are shared across candidates. This supplies no inferred measured provenance or
temperature dependence for an authored curve.

`rebind_capacitor_curves_candidate` replaces the complete curve binding set at
fixed charge and flux. Removed bindings return to their original graph C. The
candidate reports signed reversible work from the exact piecewise integral of
V_new(q)-V_old(q), summing over the union of both original knot partitions.
Voltage differences are integrated directly, so resolvable work survives below
the spacing of biased endpoint energies. The physical parameter owner pays or
receives that work when publishing the candidate; it is not resistor heat.
Ideal-source consistency, curve support, bounded work and cancellation remain
admission gates. Retained source-derived linear storage cannot claim authority
for a different authored curve. The existing quadratic gap-motion owner also
refuses these curves until a joint geometry-dependent energy law is supplied.

The ordinary circuit independently computes its stored-energy increment from
each capacitor's declared integral and the native full mutual-inductance
gradient, including represented-midpoint rounding. Compensated summation and
bounded stack scratch avoid subtracting large endpoint energies or refactoring
the inductance matrix. `MaterialCircuitStep::stored_energy_change_j` exposes the
native increment after this independent check; primary finite-enthalpy,
conductivity and piston-thermal accounts consume it. Endpoint energies remain
diagnostics. Source-free energy, independent RC trajectories across a knot,
contact/bulk heat partition, conductivity refresh, mutual-field retention,
signed parameter work, ideal-source/source-authority refusals and atomic
cancellation are exercised in `tests/material_circuit_capacitance.rs` with
synthetic material data.


### Source-backed real gas in the original inertial piston

`GasPiston<G = ClosedGasParcel>` accepts the original ideal thermochemical parcel
or the existing `IsochoricRealGasParcel` through a sealed constitutive interface.
Both use the same native wall coordinate, free/fixed support, impulse solver,
acceptance gates and port account. `PistonStep<S>` retains the complete resolved
midpoint state. Existing ideal defaults retain their original applicability.
The reciprocal finite-body heat solve, state-dependent body-power solve and
source-bound translating inertia now accept either gas owner. Specialized
reaction, double-acting and network adapters retain their stated gas
applicability. Explicit real-gas nozzle adapters retain their own source,
flow-splitting and homogeneous-phase admission contracts.

The gas owner supplies conservative heat/volume candidates, work-conjugate
midpoint pressure and its discrete endpoint-volume derivative. Actual resolved
density and sound speed supply before/midpoint/endpoint acoustic clock admission;
real-fluid mechanics never use ideal-reference density or acoustics. The same
midpoint pressure supplies force, gas work and fixed-support reaction. Accepted
gas energy, wall kinetic energy, exterior-pressure work, additional-force work
and heat must close before final cancellation admission and atomic publication
of gas, motion and the shared port window.

The PR mode retains the homogeneous frozen-composition source-admitted branch.
It supplies no phase equilibrium, gas latent heat, reacting chemistry,
friction/contact stops or exact-isentropy/nonnegative-entropy guarantee. Focused
regression sources cover dense-gas motion and energy accounting, fixed-wall
heating, real-EOS acoustic admission, and source/cancellation refusal with replay.

#### Dense-gas compression with finite wall heat and source-bound inertia

`ThermalPistonCandidate<S, G>` and `PoweredThermalCandidate<I, P, G>` retain the
original `GasPiston<G>` endpoint. Their existing bounded feedback solve evaluates
the gas at its actual conservative energy/volume midpoint and the thermal body
at its enthalpy midpoint. The same signed interface heat enters the gas and leaves
the finite body; it cancels from the independently checked combined gas, body,
kinetic and external-work balance. A body enthalpy plateau may supply latent heat
without claiming gas phase equilibrium. State-dependent body power receives those
same two material states and is counted once as external energy.

`MaterialConductionLayers` consumes the existing `FluidTransport` snapshot, so
real-fluid density, heat capacity and conductivity stay attached to the resolved
dense-gas state. Body conductivity retains the original card, selection policy,
auxiliary axes, temperature and optional liquid-fraction coordinates, with typed
source refusals. No ideal-reference transport or name-based property fallback is
introduced.

`MaterialGasPiston<G>` obtains translating mass from the original `ResolvedBody`
inertia and retains that body together with the original native piston. Its
thermal and powered candidates preserve both owners and the native port account;
they do not reconstruct a fresh piston or silently change mechanical mass when a
thermal body heats or changes phase. Source-bound layer conduction and actual-
midpoint power therefore continue accepted compression without dropping density
receipts, resetting motion, or counting source energy twice. Thermal storage may
represent a separate wall, and its uniform-temperature, fixed-geometry and
negligible-deformation-work assumptions remain explicit.

`MaterialGasPiston::thermal_layers_candidate` connects those original conduction
layers to the retained material-bound mechanism for either admitted gas owner.
Its `PortedModule` implementation forwards every schema, coordinate chart,
junction, accepted window, transfer, observation error and audit to the original
native piston. It introduces no second port clock or reconstructed exchange.

Focused regression sources cover repeated dense-gas compression with sourced
finite-wall heat, latent-wall cooling, actual-midpoint state-dependent power,
source-derived inertia, fixed supports, isolation, and refusal/cancellation
without accepted-state or port-history mutation. These extend the existing
integrators and acceptance budgets; they do not introduce a higher-order scheme
or an unconditional stability guarantee.


## Current tetrahedral-surface contact geometry

`finite_solid::contact::surface::FiniteSolidContactSurface` prepares the original
`ClosedTetPressureSurface` exterior, binds its exact reference coordinates and
connectivity, retains canonical face/tet incidence, and caches the query owner's
triangle BVH. Queries refit current coordinates without rediscovering topology.
They retain original tetrahedral barycentric material points and actual rigid
sphere arms, including offset sphere centres and rigid orientation. Shared edge
and vertex witnesses deduplicate by stable topology-derived contact IDs;
permuting sphere inputs does not change the canonical result.

The actual sphere and solid points remain distinct at positive gap. No common
point is fabricated to suppress the angular couple from a tangential impulse.
The complete geometry feeds the existing consistent-mass/full-inertia
`rigid_contact_candidate`, preserving off-centre torque, physical frictional
spin, and its original spin/orbital angular-momentum and work audits.

Raw `sphere_separations` receipts retain negative gaps for event bracketing.
`sphere_contacts` and `admitted_contacts` refuse penetration without changing
either configuration. Inverted or collapsed current tets, interior or unresolved
sphere centres, reference mismatch, invalid/range-limited geometry, duplicate
sphere IDs, identifier overflow, cancellation and exhausted work/output bounds
refuse. A sphere centre outside the refitted whole-surface AABB skips the
per-tet containment scan; other centres receive the complete bounded test.

This route describes discrete P1 geometry. The original volume mesh must be
conforming and globally non-overlapping; local tet admission does not certify
global injectivity. No arbitrary-flow CCD, penetration repair or physical
impact-loss heat partition is inferred from the geometry query.

`tests/finite_solid/contact/surface.rs` covers the deformed material trace,
separate positive-gap witnesses, canonical edge/vertex deduplication, an actual
off-centre full-inertia impact, Coulomb-induced spin of a homogeneous sphere,
and penetration/interior/budget/cancellation refusals. The internal current-tet
fixture checks inversion and collapse. These are authored numerical oracles,
not experimental material validation.

## Native dielectric polarization and finite thermal feedback

`DebyeCircuit::rebind_source_candidate` refreshes the complete original dielectric
source at fixed individual instantaneous/polarization charges and original flux
linkages. It admits the same material claims, pole roles, nonthermal context and
fixed geometry. Stable fixed-charge free-energy change is returned as signed
parameter work. The optional complete conductor binding retains original branch
indices, conductor source identity and geometry across subsequent R(T) updates.
Internal polarization friction remains a dielectric process. Incompatible ideal
voltage clamps refuse instead of projecting charge or inventing impulse work.

`circuit_enthalpy::dielectric::step_candidate` couples that native history to the
existing finite caloric, conduction and Newton--Krylov owners. Every trial starts
from the same accepted charge and thermal histories, resolves all original
material resistors at the trial endpoint, and evaluates the original Debye source
at the dielectric cell's actual temperature. It uses the existing native circuit
midpoint step; there is no separate dielectric integrator. Per-branch native
resistor heat goes to each original conductor's receiver, while the sum of native
polarization losses goes once to the dielectric's own finite caloric cell.

The material law supplies Helmholtz free energy psi, entropy S, internal energy
U=psi+T*S and a conservative fixed-charge heat-capacity bound. The distinct signed
reversible transfer is Q_rev=-delta(psi_parameter)-delta(T*S). It enters the
original caloric residual separately from irreversible polarization/resistor
heat. The final physical account removes this internal transfer from external
heat and checks the complete electrical internal-energy plus caloric increment
against electrical supply and actual external heat. Parameter work and delta(T*S)
use their stable material increments; ordinary endpoint differences do not own
the field transfer. Source receipts and both candidate histories publish together.

One complete uniform dielectric volume must coincide with one fully solid,
constant-density caloric cell from the same original material card and context.
The admitted temperature interval needs full caloric support. Explicit bounds
on every charge admit strictly positive combined caloric-plus-field capacity
throughout that interval, and are rechecked after each electrical trial.
Temperature-dependent relaxation times may vary along supported positive source
paths; capacitance paths use the material thermal law's affine interval. Other
original capacitor/inductor coefficients remain fixed. Multiple receivers and
native thermal links are supported, but this owner represents one dielectric
branch and one uniform dielectric temperature. Phase changes, geometry motion,
ideal-source impulses, nonlinear polarization saturation and breakdown are not
supplied. Endpoint material coefficients make the general coupled scheme first
order; energy closure is not an unconditional entropy or temporal-error proof.

`tests/dielectric_thermal.rs` independently checks midpoint polarization and total
internal energy for both signs of C(T), absorption self-heating, the zero-pole
external-heat limit, tau(T) temporal refinement against a scalar reference,
distinct finite dielectric/conductor receivers with R(T), and atomic capacity,
work-budget and cancellation refusal. `tests/dielectric_rebind.rs` exercises
conserved native histories, signed parameter work and source/ideal-clamp refusal.
All property data in these tests are explicitly authored numerical fixtures.

`contact::step::free::rigid_free_flight_candidate` advances the original finite
solid objective energy/momentum discrete-gradient owner and each original rigid
`step_verlet_with` owner over one positive interval. It retains native state receipts and returns both candidates
only after complete momentum and work checks. Solid forces come from the same
prepared external-load owner as the native solve. Rigid gravity is included
once in force work; a solid may bind its own explicit stationary field, whose
native nodal loads, potential work and integration defect retain separate accounts.

The solid account uses consistent `M delta-v`, current midpoint coordinates and
the full angular increment `x_average cross delta-p + delta-x cross p_average`.
Original stationary-support reactions contribute external linear and angular
impulse and zero work. The rigid account uses actual endpoint positions,
attitudes and momenta, both endpoint force moments and rotated body-torque
half-impulses. Native half-kick momenta and full inertia determine torque work;
quaternion polarization resolves the actual spin change without subtracting
two large world spins or assuming free-drift conservation.

Stable kinetic changes, original strain changes, source work and signed solid
and rigid-drift integration defects remain distinct. The absolute sum of the
individual numerical defects is bounded, so opposite defects cannot hide each
other. Momentum and work gates also bound the sums of individual-body residual
magnitudes, so opposite lost kicks cannot cancel in the complete-system account.
Unrepresentable kicks refuse. Inputs remain unchanged on native failure,
cancellation or late audit rejection. This entry does not admit a collision
path, certify trajectory accuracy, model exact externally forced Euler flow,
or infer heat. The no-pressure entry retains the native refusal for a
pressure-bound solid.

`rigid_free_flight_candidate_with_pressure` holds an explicit uniform pressure
on the original closed material boundary. The original pressure-aware discrete
gradient advances material energy plus pressure-volume potential. Its native
`-p Delta V` work stays visible in `report.solid.pressure_work_j` and enters total
source work exactly once. Closed uniform pressure adds no net external force or
torque; its algorithmic fixed-support reactions remain in the original support
impulse account. Nonuniform face pressure is not supplied by this entry.

Focused G1/G3/G4 fixtures exercise constant acceleration and axial spin/work,
deforming-solid and asymmetric-rigid conservation, stationary support reaction
moments, and lost-kick/cancellation refusal. Material coefficients are synthetic.

### Automatic transient rigid-sphere/finite-solid trajectories

`finite_solid::contact::step::sphere_trajectory_candidate` composes the original
current-surface detector, simultaneous consistent-mass/full-inertia contact
solver, objective energy/momentum solid motion and native full-inertia rigid
motion. Each rigid body supplies one COM-centred sphere. The entire trajectory remains a private candidate until all requested
time, geometry, physical balance and work gates succeed. Current contacts are
resolved together at their actual positive activation gaps; original material
witnesses are never moved to eliminate tangential angular couples.

The prepared triangle BVH refits unions of native endpoint coordinates. Directed
separating-plane bounds admit every declared straight sphere/triangle path;
interval tetrahedron determinants exclude intermediate element collapse.
Rigid-sphere pairs use the same directed separation authority with an upward
radius-sum enclosure. A possible rigid/rigid encounter refuses because those
rows are not yet part of the shared rigid/solid complementarity problem. No
pairwise impulse substitute is applied. An unresolved solid-contact window
retries a shorter original native step from the same unpublished post-impact
state. Failed candidates do not advance any caller-owned q/v/a, attitude,
momentum, time or accepted history.

The declared straight coordinate interpolation is the numerical collision path,
not an enclosure of the true nonlinear flow. Positive activation distance,
maximum timestep and work caps are explicit refinement controls. Offset sphere
centres need a curved rotational path and remain limited to discrete event
queries. Sustained loaded contact, unresolved grazing and Zeno accumulation can
exhaust progress or work bounds; they never create artificial rest or runout.
The original no-pressure entry retains the explicit-pressure refusal. This is
not a general Moreau--Jean trajectory, global mesh-injectivity proof or a
claimed global convergence order for impacts.

`sphere_trajectory_candidate_with_pressure` and
`mixed_sphere_trajectory_candidate_with_pressure` admit an explicit uniform
pressure on the original closed surface before contact work. Both reuse the
same event/retry core and preserve that held pressure through every shortened
trial. Each accepted free phase uses the original pressure-aware objective
solid step. Impacts retain the original positions, so fixed-pressure work
telescopes across the whole call. Simultaneous contact loss, source/pressure
work, support impulse and signed/absolute numerical defects stay separate.
These entries do not extend the sustained loaded-contact lane or infer a
nonuniform face-pressure law.

Whole-call accounting retains stable kinetic/strain increments, original source
work, stationary-support impulse and moment, signed/absolute smooth integration
defects and total unallocated contact loss. Absolute local balance residuals
accumulate before global cancellation. Native phase and endpoint clocks must
agree with requested elapsed time within the declared binary64 rounding gate;
unrepresentable timing refuses instead of overwriting native state. Gravity
work is included once: the rigid integrator owns its gravity, while the solid's
bound gravity owner retains its own sampled loads and potential/work account. No loss becomes heat
without an independent physical partition.

Focused G1/G3/G4 fixtures specify an automatically detected elastic bounce,
friction-driven spin and support moment through real attitude evolution,
otherwise unseen sphere-pair crossing, late budget/offset/cancellation rollback,
and loss of requested time at a large native timestamp. These are authored
analytic references, not experimental validation.

The default-DG/pressure extension adds four public regression cases for nonlinear
free flight, a freely recoiling pressurized impact, pressure-preserving geometry
retries/zero-pressure parity and binding/cancellation refusal. The full Cx-bound
public target remains uncompiled and unexecuted where the runtime dependency
closure is absent. A bounded actual-native numerical fixture executes original
consistent-mass assembly/FGMRES initialization, simultaneous contact, material/
pressure stress and tangent, time integration and rigid impulses/drift. Its five
finite-step comparisons reduce a free-spin work defect from about 7.68e-5 J to
1.00e-15 J and a post-impact defect from about 2.53e-7 J to 1.89e-15 J. The
pressurized impact has 0.375 N s impulse, 0.09375 J physical contact loss and
0.0156686896022 J pressure work, closing work within 3.47e-15 J. This direct
composition does not execute the public geometry/retry or receipt wrappers.
The existing native frictional-impact regression executable also passes 18/18;
these checks are not a full Cargo/DSR gate or experimental qualification.

### Ordinary heavy solids on the resolved immersed-fluid owner

`material_immersed::settling::{SubmergedMaterialSpheres, SubmergedMaterialEllipsoids}`
bind original ordinary solid sources to the same native D3Q19 moving-body and
sphere-film solver. The analytic shape supplies displaced reference-fluid mass
and the complete rotated tensor. The carrier receives their explicit difference
from the original physical mass/inertia; source cards and receipts stay unchanged.
One frozen fluid snapshot supplies density, units and bulk/film viscosity.

The grid uses reduced pressure with zero fluid gravity. Each carrier receives
`(M_s-rho_f V)g` once, plus separately authored center loads. Strictly positive
excess mass and full physical excess tensor are required; equal/lighter solids
refuse without an artificial mass. Construction releases all bodies from rest,
matching the quiescent fictitious interior fluid. Full immersion, fixed source
geometry, native clearance/Mach bounds and existing film assumptions apply.

The native solver remains the sole time/pose/population owner. Its raw receipt
retains represented-fluid/excess-carrier impulse/work and numerical rotation
error. Optional SI observations reconstruct dry-body momentum using displaced
inertia acceleration and Archimedes terms; these are distinguished from measured
population reactions. No fallible conversion follows native publication.
Density-preserving source viscosity refresh retains accepted history and updates
bulk flow and any films together.

Optional current-state diagnostics reuse native sampler scratch to report
occupancy-volume error, actual interior velocity slip and momentum mismatch.
They neither enforce interior rigidity nor alter state. Continuum excess-inertia
equivalence requires rigid interior fluid; the finite partial-saturation law and
subcell quadrature approximate it and require spatial/time refinement. This is
not a free-surface, arbitrary-density, calibrated-drag or thermal-dissipation
closure. No complete fluid-energy theorem is inferred from body work.

Focused synthetic tests independently check source mass/tensor subtraction,
actual population-plus-excess momentum under gravity, reconstructed dry momentum,
represented work, rotated ellipsoids, simultaneous sphere-film reciprocity,
density-preserving source refresh and atomic native/source/cancellation refusals.

## Qualified causal evanescent reflection

`mm_line::evanescent::MatchedEvanescentPorts` now realizes the independently
matched higher radial ports of a lossless circular guide through physical LC
storage. For each native radial cutoff `wc`, its source is
`Z(s) = z0*s/sqrt(s*s + wc*wc)`. A positive spectral quadrature gives a series
chain of parallel inertance/compliance cells with `L=z0/(N*wc)` and
`C_j=N*cos(theta_j)^2/(z0*wc)`, `theta_j=(j+1/2)*pi/(2*N)`. Static inertance is
exact. The original gas, radius, roots, ordered cells and resolution remain
inspectable; no pole fitting, finite-resistor substitute or gain repair occurs.

The existing `fs-vfit::relaxation::ParallelCellLoad` owns both physical states
and midpoint stepping. The instantaneous wave closure solves
`(D+z0)*U = 2*a-history` on the current sample and returns `b=a-z0*U`.
Its discrete identity is `delta H = dt*(a*a-b*b)/z0` for any representable input
sequence, in exact arithmetic. Reactive storage can return energy without
amplitude clipping. All higher modes prepare before any state commits; a late
mode refusal preserves every accepted state and observation. At most seven
ports of sixteen cells are retained; sample stepping allocates nothing. The
affine next-sample response is available to an implicit source junction.

Physical accuracy has an authored band, not an all-frequency claim. Its upper
frequency must remain below Nyquist and, after the midpoint/Tustin warp, below
every retained cutoff. For normalized warped frequency `r`, the relative
impedance quadrature error is `e=2*q^(2*N)/(1+q^(2*N))`,
`q=(r/(1+sqrt(1-r*r)))^2`. The complex reflection quadrature bound is
`2*r*e/sqrt(1-r*r)`. The additional time error is the chord length between
`-exp(-2i*asin(omega/wc))` and its warped-frequency value. Their sum bounds
complex reflection error throughout the band and must meet the caller's
budget. These analytic approximation bounds exclude floating-point roundoff.
The discrete LC model remains passive outside the band, but it does not
reproduce above-cutoff propagation/radiation or viscothermal loss.

`UniformModalLine` supplies the useful bounded composition: one positive-length
lossless cylinder, the actual original one-mode `MmLineBank` FIR owner, and
these independently matched higher ports. The plane response keeps its original
delay, arithmetic and finite-grid admission receipt. Scalar radiation and strict
tabulated loads keep their disclosed matched-higher-mode closure. Nonuniform
geometry, Bessel loss, full flanged modal coupling, closed/pressure-release ends,
matched plane loads and diagnostic fallback sources refuse explicitly. There is
no silent replacement in `MatrixFirLine`, no extra processing delay, and no
claim that this bounded rung completes a general multimodal radiation realizer.
Regressions compare the native modal authority, separate three-rung quadrature
and clock ladders, sampled harmonic response, independently reconstructed LC
storage, wave work, same-sample closure, rollback, and the previously refused
uniform resistance fixture with bit-identical scalar controls.

## Simultaneous creep and plastic heat on spatial enthalpy

`solid_enthalpy::creep::{CreepCell, CreepSolid, CreepCandidate}`
accept the material owner's `CreepSourceFamily`, defaulting to their original
Norton source. Selecting `NortonJ2Source` evolves both native creep and plastic
history through the existing tetrahedral solve and the same thermal feedback
loop. Actual source-selected yield stress may change with local temperature;
E, nu, energetic H and reference density remain fixed over the admitted span.
Both physical dissipation channels enter the caloric source exactly once;
hardening storage and numerical loss remain separate.

Original caloric identity, reference mass/geometry, conductances, reservoirs,
observer hooks, actual-temperature re-equilibration and atomic publication are
shared. `tests/solid_creep_plastic_enthalpy.rs` crosses a decreasing sourced
yield curve under load, checks actual permanent deformation against an
independent scalar thermal/mechanical reference, and exercises nonconvergence,
observer failure, cancellation and retry with both memories intact. The
unwrapped source selects zero expansion. Variable stored-energy coefficients,
finite deformation, phase flow, damage and life qualification remain outside
this reduction.

## Reversible nonlinear dielectric gap motion

`material_circuit::motion::dielectric::ReversibleDielectricGapMotion` advances a
source-backed nonlinear dielectric branch, the original circuit and a guided
material body in one native descriptor. The retained reversible material owns
the monotone D(E) source and its inverse and integral. With fixed electrode area
A, gap d and conserved free charge q, the joint potential is
H_field=A*d*w(q/A), w(D)=integral E(D) dD. Terminal voltage is d*E(q/A); the force
on increasing gap is -A*w(q/A), the derivative at fixed charge. Both charge
polarities use this same even energy and attractive force. The nonlinear energy
is not inferred from a secant or differential capacitance or from half qV.

Body mass comes from the original regional density/geometry compiler. An
optional original `AxialSupport` contributes its material restoring force,
elastic potential and consistent modal inertia. Gap and momentum extend the
original charge/flux coordinates with a canonical power port; the original
solver, ideal-source/transformer topology, bulk conductor and series-contact
loss owners remain in use. Other authored capacitor curves and the complete
signed mutual-inductance matrix retain their original native energy. Immutable
prepared dielectric and capacitor source tables are shared across proposals.
The `routing_circuit()` accessor is explicitly a topology/source view: the
selected graph C is only a routing placeholder, and does not represent the
active field's physical storage. Advancing or auditing the isolated routing
view does not substitute for the joint owner.

Independent electrical, kinetic and support energy increments use the original
constitutive integrals and represented-midpoint corrections with compensated
summation. The combined native storage supplies the corresponding analytic
midpoint remainder, including the charge-to-displacement transform; it does
not subtract large endpoint energies. Absolute endpoint energies remain
diagnostics. Each candidate checks these increments against native work and
dissipation, retains separate original bulk/contact heat, and reports external
force work from force times actual gap displacement. Field work is one internal
exchange between the electrical and mechanical stores, never additional heat.

This bounded model has one uniform, fixed-area moving dielectric gap. Its
source must depend only on field magnitude; an actual temperature-dependent
source slice refuses because a frozen free-energy law does not supply missing
reversible heat. Applicability-only temperature bounds are retained. A live
thermal support likewise requires its separate reciprocal thermal owner and
is refused here. Original field/charge support, open positive gap limits,
support strain limits, finite controls, bounded descriptor capacity and
cancellation are admission gates. Native initial and endpoint consistency must
preserve every charge, flux, gap and momentum coordinate; incompatible ideal
clamps refuse without a charge projection or unreported impulse. Both circuit
and body histories publish together only after the complete physical audit.
No fringing, electrostriction, material-flow work, changing area, dielectric
memory, breakdown, collision law or thermal evolution is supplied by this
uniform reversible profile. Energy closure alone is not a temporal-error proof.

`tests/dielectric_motion.rs` exercises independent fixed-charge energy
derivatives and constant acceleration for both polarities and source masses;
nonlinear RC motion across a material knot against independent RK4 and step
refinement; the existing linear dielectric-motion limit; retained unrelated
nonlinear capacitance and full mutual field storage; original contact heat;
source support inertia and elasticity; tiny SI work below a large unchanged
capacitor bias; and source/domain/ideal-clamp/cancellation refusal. The source
curves are authored numerical fixtures, not claims of measured material data.


## Source-selected thermal Voce response in the loaded solid

`ThermoplasticSolid::with_voce_hardening` accepts one optional independently sourced saturation selection per tetrahedron before mechanical evolution. Each selected element retains the original material/caloric card, pinned Q(T)/b claims and complete source bundle. `None` leaves that element unchanged; unsupported simultaneous Voce/overstress selection refuses through the material owner.

The existing deposition heat timeline, spatial Newton solve and accepted J2 histories remain the only evolution owners. Actual endpoint temperature selects the native Voce stress and tangent. The existing mechanical audit includes source-coefficient storage, saturation entropy and signed reversible exchange, while physical yield heat and backward-Euler numerical loss remain separate. Failed or unaccepted candidates do not alter caller-owned thermal or mechanical state.

The heated-member regression uses six connected tetrahedra with held axial grips and three physical heating intervals. An independent scalar axial equilibrium and thermal balance oracle verifies the temperature, all element stresses/history, reaction force, lateral plastic contraction, original source receipts and fixed reference mass. Affine Q(T) softening lowers held-grip tension and increases retained plastic strain relative to the constant-Q control. Separate formulas verify mechanical work, free-energy change, coefficient storage, delta(T s), physical dissipation, reversible caloric heat and numerical loss. Per-element selection and unsupported/late/missing binding refusals are covered. These synthetic fixtures verify the coupled implementation; they do not calibrate a physical material or add thermal Chaboche, finite-strain or hexahedral heat transport.


### Source-backed expansion and reciprocal heat in spatial permanent creep

`solid_enthalpy::creep::CreepSolid` is the shared owner for the original
zero-expansion Norton/J2 sources and their `ExpandingCreepSource` composition.
This replaces the former `ZeroExpansionCreepSolid` API name; all in-repository
callers use the common owner. It retains the existing native tetrahedral solve,
actual caloric-cell mapping, conduction/reservoir/event owner and one atomic
accepted mechanical/thermal/observer state.

Construction verifies the original source/caloric identities and fixed reference
mass. For an expanding source, the zero-displacement virgin birth must be at an
actual stress-free temperature. Every piecewise-linear segment of the actual
enthalpy chart, with downward-rounded capacity arithmetic, must keep the sum of
caloric capacity and the source's full-interval mechanical capacity floor positive.
A valid coefficient set with insufficient combined capacity is refused.

The original feedback loop now transfers physical creep/plastic heat plus signed
reciprocal thermal exchange. Accepted history retains the old cell temperature;
all trial solves restart that same state, and final mechanics are resolved at the
actual returned temperatures. Per-cell absolute heat defects, temperature gates,
thermal energy and the combined free-energy + T*entropy + caloric balance all
precede publication. This supports cooling during stress relaxation while creep
dissipation stays positive. Numerical loss is excluded from heat.

The existing `solid_creep_plastic_enthalpy` target retains its prior joint-law
cases and adds stress-free heating/expansion and restrained heating that yields
and creeps before relaxing with negative reciprocal/net caloric transfer. The
references solve independent scalar force and heat balances and verify rollback
from nonzero histories after observer failure. The model remains quasistatic,
small-strain, fully solid and fixed-reference-geometry, with constant E/nu/H and
reference density. Finite deformation, moving-mesh transport, variable elastic
storage and service-life qualification are not implied.

The coherent native check refreshed 15 affected production crates and typechecked
all six frame/tetrahedral/coupling integration targets, including both
`solid_creep_enthalpy` and `solid_creep_plastic_enthalpy`. Native integration tests
were not numerically executed. The separate 56-test material behavior suite ran
successfully with genuine dependencies; neither check is a full Cargo/DSR gate.

### Original-source immersed solids at light, neutral and heavy density

`material_immersed::dry::ResolvedImmersedBody` binds one original spherical or
rotated ellipsoidal solid to the native implicit dry-body/interior-fluid owner.
The existing material/geometry validation, source receipts, SI/lattice scaling,
complete positive dry inertia and viscosity-rebind path are reused. Analytic
original volume supplies reduced gravity `(M_s-rho_f V)g` once; sampled occupancy
is never substituted into that load or renormalized. Initial body momentum is
retained even at equal density. Source/fluid density changes through rebind remain
refused because they require conservative transfer.

A successful step returns the unchanged native physical/momentum diagnostics in
lattice units. Unit conversions and reduced-gravity potential are separate
observations; no fallible post-publication conversion can report a rollback that
did not happen. Rejected implicit steps preserve fluid, dry body, original source
and accepted-step count. One fully immersed body only: stationary/periodic fluid
boundaries, no contact film, no free surface, and no exact rigidity/GCL or
unconditional light-body stability claim.

### Original-source physical dry assemblies and sphere films

`material_immersed::dry::ResolvedImmersedBodies` extends original material
ownership to one through eight simultaneous light, neutral or heavy bodies.
Its sphere and rotated-ellipsoid constructors retain every source, analytic
original volume, complete positive dry inertia and authored initial motion.
All bodies share one source-owned reference fluid and lattice/SI scaling.
Reduced gravity is `(M_dry-rho_reference*V_original)*g` once per body;
the reduced-pressure fluid has no additional gravity source.

Sphere assemblies use the original native implicit dry-body and simultaneous
lubrication owner. Bulk and normal/optional tangent-spin film resistance obtain
viscosity from the same fluid snapshot. Film loads participate in every coupled
candidate, including the measured streamed interior-fluid momentum correction;
they are not applied after the solve. Rotated ellipsoids use the same original
full-tensor dry-assembly owner with native shape axes and enclosing-sphere
clearance. No ellipsoid film law is inferred.

`step_forces_si` accepts world center forces and couples and returns the
unchanged complete `DryBodiesExchange3` in lattice units. Explicit observation
conversions preserve the native separation of external, resolved-fluid,
interior-correction and film accounts. Density-preserving rebind updates bulk
and film viscosity together without discarding accepted fluid/body history.
Density changes remain refused because they require conservative transfer.
All unit conversions and the next step count are checked before the native
atomic publication; refused calls preserve accepted sources, populations,
bodies and time. Prepared sampler/solver scratch may change on a refusal.

The bounded applicability remains stationary or periodic fluid boundaries,
fully immersed nonoverlapping shapes, native clearance/Mach/volume-change
limits and a converged finite dry-body correction. Positive dry inertia alone
does not guarantee convergence at every density ratio or timestep. Exact
rigidity/GCL, free surfaces, deforming solids, a thermal destination for film
loss and unconditional light-body stability are not supplied.

Four public regression groups cover mixed source densities and original inertia,
independently rebuilt full-fluid/dry-body momentum at a fixed origin, physical
normal-film resistance and external wrench conversion, source-preserving rebind,
atomic refusal/exact retry and rotated full-tensor ellipsoids. They are authored
and statically reviewed, with pinned formatting checks, but remain uncompiled
and unexecuted: the narrow Cargo command fails loading the missing asupersync
sibling before Rust compilation. No public-target pass is claimed.

### Material-owned dry assemblies with prescribed film surfaces

`ResolvedImmersedBodies::with_prescribed_planes` attaches SI
`PlaneFilmSi` inputs to a spherical material assembly and returns
`ResolvedImmersedPlanes`. Original material identity, density, complete dry
inertia, populations, accepted pose/momentum and step count survive attachment.
Plane normals are normalized by the original native owner. Wall-normal and
optional tangent/parallel-spin cutoff/roughness scales remain explicit and
independent of sphere-pair scales. Viscosity comes only from the same resolved
fluid snapshot used for bulk flow and every film.

One original implicit dry-body solve includes all sphere-pair and plane
resistance. The complete `DryLubricatedPlanesExchange3` remains in lattice
units; each plane separately receives its linear impulse, angular impulse
including reaction couple, and work. Its prescribed reservoir supplies the
opposite to the represented fluid/body system. Resolved-wall exchange,
authored loads, reduced weight, fictitious-fluid correction and modeled film
loss remain distinct. No reaction is applied to populations a second time.

The wrapper owns accepted plane positions. Their native constant translation
advances exactly once per accepted interval; all next-plane positions, body
load conversions and the next count are checked before native publication.
Rejected probes, cancellation and exhausted budgets retain all accepted
plane/fluid/body history. The reusable next-plane buffer is prepared at
attachment; density-preserving viscosity rebind retains geometry and updates
bulk, pair and plane resistance together. Read-only assembly access prevents
advancing body time independently of its attached surfaces. Empty planes
retain the original sphere-only result exactly.

These surfaces supply missing-film resistance only. They do not create or move
a resolved no-slip boundary, and the scenario must supply a matching resolved
boundary separately. Finite plane mass/recoil, rotational motion, contact,
ellipsoid films and thermal allocation of film loss remain outside this owner.
Original dry solve, low-Mach, positive-gap, sweep, film/mode and at-most-64-plane
admission remain authoritative.

Three public regression groups cover independent exterior-fluid/dry-body/plane
momentum balance, normal-film SI resistance and reservoir work, exactly one
geometry advance, empty-plane parity after existing history, and viscosity
rebind/refusal/exact retry after native probes. They are authored and
statically reviewed; pinned formatting passes. The public target remains
uncompiled and unexecuted because the constellation preflight finds missing
siblings and Cargo fails on asupersync before Rust compilation.

### Material-owned moving finite film surfaces

`ResolvedImmersedBodies::with_moving_surfaces` binds mixed
`MaterialFilmSurfaceSi::{Prescribed, Finite}` descriptors to the original dry
sphere/fluid assembly. `ResolvedImmersedSurfaces` retains the source descriptors,
material identity, complete source inertia, fluid snapshot, accepted state and
clock. Finite carriers use the same shared SI/lattice full-inertia conversion
as sampled bodies; body-frame points and normals move with the original native
rigid pose. Prescribed descriptors include the existing axial-spin law.

One original native dry-body correction and simultaneous resistance solve
includes every sphere and surface. Finite recoil influences all other film
responses and the streamed fluid candidate during each nonlinear probe.
Accepted surface state advances once with the spheres and fluid; failed
probes and refused steps preserve every accepted state and the step count.
All SI load conversions and the next clock value are checked before native
publication. Complete authored surface loads are world center forces/couples;
prescribed slots require exact zeros and receive no extra imposed impulse.

Read-only observers expose accepted finite-carrier SI states, actual world
plane geometry and native local point velocity. The assembly's original mass,
impulse, angular-impulse and energy converters interpret the unchanged native
report. Prescribed actuator work, finite-carrier kinetic storage, native
free-rotation energy defect and isothermal film loss remain distinct. The
sphere/exterior-fluid account is an open subsystem; its surface reaction must
not be mistaken for a second external load on a finite carrier.

The wrapper delegates source-consistent viscosity refresh, initial physical
flow and held uniform/spatial fluid-force rebinding to the original assembly.
No independent surface viscosity or mutable assembly escape is introduced.
Density changes keep the original conservative-transfer refusal. Empty and
prescribed-only inputs retain the original admitted physical trajectories.

Each finite descriptor is an independent local planar missing-film carrier.
Source inertia does not certify an infinite plane or its geometric extent.
The caller owns matching resolved-wall geometry and the applicability of the
planar approximation. No sampled carrier interior, buoyancy, inferred weight,
shared multi-face carrier, contact, exact swept collision detection, moving-wall
GCL or thermal destination is supplied. Complete finite-carrier external loads
must therefore be authored explicitly; native admission and dry budgets apply.

Four public regression groups cover source-derived full inertia and SI/local
velocity, repeated mixed-system momentum/work and accepted motion, exact legacy
continuation, and live coefficient refresh with refusal/cancellation/retry.
They were authored but remain uncompiled and unexecuted: runtime compilation
was killed before the changed native crate, and the execution environment
subsequently became unavailable.

### Physical initial flow and external fluid loads on material bodies

`RetainedFluidBodies`, `ResolvedImmersedBody`,
`ResolvedImmersedBodies` and `ResolvedImmersedPlanes` expose the original
native uniform-flow and held-force owners in SI units. The shared source-scaled
container performs all conversions; the dry/plane wrappers retain their existing
body/film solver, original dry inertia, material identity and accepted clock.

`initialize_uniform_flow_si` prepares reference-density physical fluid velocity
only before the first accepted interval. It leaves authored body momentum and
prescribed plane motion intact and checks the source-scaled scenario speed bound.
The original equilibrium is transformed into the CURRENT Guo forcing frame.
A body/wall velocity mismatch is a startup transient, not a solved steady flow
or exact rigid interior. Initial fluid preparation is not an evolved-state
impulse and is refused once time has advanced.

`rebind_uniform_fluid_force_si` accepts an external fluid-only volume force in N/m^3.
Conversion is
`F_lattice = f_SI dt^2/(rho_reference dx)`, using the same retained source
density as body/fluid momentum units. The native uniform-force reframe preserves
physical mass, momentum and force-corrected stress without advancing time or
supplying work. Subsequent original steps apply the held force. Dry-body
interior-fluid correction already removes the fictitious interior source;
callers must not repeat it in solid loads or add the removed reference gravity.
The original native receipt owns external-fluid impulse/torque and its declared
local work account. Viscosity rebind retains the held force and mechanical state.
An affine pressure `p=p_periodic-f.x` additionally requires the separately
authored dry-body traction `V_original*f`; periodic perturbation pressure cannot
recover it. That physical solid force is distinct from the fictitious-fluid
subtraction and is not a duplicate fluid deposit.

`rebind_fluid_force_field_si` accepts a dense x-fastest field in N/m^3,
one vector per original grid cell, and adds it to the retained uniform load.
It selects a piecewise-held spatial history on the same native owner. All SI
conversions are admitted before native publication; static solids receive no
samples, while partially occupied body cells retain the field. The existing
dry solver subtracts their midpoint interior force and torque exactly once.
Rebinding preserves physical density, momentum, full force-corrected stress,
body/plane history and accepted time. A zero field disables only the spatial
component. The spatial proposal allocates one native force field and reuses
collision scratch without a second converted field. Dry single-body, assembly
and prescribed-plane owners forward to this same operation.

There is no mutable grid escape hatch,
pressure-inlet/outlet closure, variable-density law or finite driving reservoir.
The older excess-inertia settling approximation retains its quiescent-release
scope; these arbitrary-flow ports are not added to that approximation without
a corresponding interior-momentum derivation.

Three public regression groups cover independent SI fluid-plus-additional-body
momentum under an applied load, initialization at nonzero held forcing,
momentum-preserving source switches, source-viscosity updates, dry/spherical/plane
forwarding and initial-time/invalid/cancellation refusal. Pinned formatting and
static review pass. They remain uncompiled and unexecuted as a complete public
target: the local constellation preflight finds missing siblings and Cargo
refuses the absent asupersync manifest before Rust compilation.

### Continuous modal reflection on the complete exact flow

`mm_line::continuous::ContinuousModalLine` is a full-state lossless
smooth-profile reflection runtime. It consumes the original
`fs_duct::modal::time_domain` source, retains every canonical state, and
connects every inlet mode to ONE `PreparedMultiportScatteringZoh` instance.
For original modal pressure P and inward volume flow U, every real wave
reference is Z0=rho*c/S_in: U=(2a-P)/Z0 and b=P-a. This reference is a wave
coordinate choice, not a falsely matched higher-mode characteristic. The
physical coupling, transverse stiffness, propagating and evanescent response
remain in the source. No scalar closure, fitted damping, removed state,
extra processing delay or midpoint/Cayley propagator enters this path.

`reflection_relation_into` exposes the current-interval mean reflected wave
as a complete affine N-port relation, including off-diagonal entries.
`push` advances the whole jointly closed system through the shared exact
held-wave propagator and commits history only after every candidate state,
output and independent port-energy account succeeds. Warm calls reuse
prepared memory. Reflected output is the interval MEAN b, not a pointwise
waveform or RMS amplitude. Original per-port integrals of b(t)^2/Z0 remain
separately available; replacing them by dt*mean(b)^2 loses within-interval
variation. Reference exchange is boundary work, never material absorption.

The spatial image has rigid or pressure-release far ends, lossless fixed
coherent gas and the original certified m=0 basis. Spatial and radial
refinement remain independent of exact linear propagation; input-hold error
for a varying source still needs a clock ladder. This bounded rung does not
complete general modal radiation, nonlinear-junction qualification, bead .88,
or an all-profile frequency-error certificate. Existing full-state exact-flow
capacities remain unchanged; requests that exceed them refuse.

`tests/continuous_modal_line.rs` drives coherent nonzero waves at every
retained inlet port, checks same-window cross-mode reflection and the
original energy balance with zero fictitious heat, and checks each integrated
outgoing energy against its mean-wave lower bound. Exact held waves are
invariant across 1/2/4 interval subdivisions. A separate 32/64/128 clock ladder
refines a piecewise-held incident waveform obtained by midpoint SAMPLING of a
prescribed smooth input; system propagation remains exact ZOH on EVERY
interval. Wrong shapes, nonfinite input and finite-overflow candidates leave
all accepted state, means and port accounts unchanged.

## Source-backed stationary nonlinear capacitors

`material_circuit::capacitance::dielectric::MaterialDielectricCapacitor` retains
one complete `ReversibleDielectricMaterial`, original branch and fixed physical
area/gap. Every original D(E) knot maps directly to `q=A*D` and `V=g*E`; the
zero-remanence isotropic law supplies both polarities. The original native
`CapacitorChargeCurve` owns its exact piecewise integral, voltage, discrete
gradient and charge history. No fitting, resampling, second store or alternate
integrator is introduced. The native signed-knot budget is enforced before
allocation. A source with actual external-coordinate dependence refuses: a
frozen thermal free-energy curve alone cannot power a closed ordinary circuit.
Applicability-only conditions and their original source receipts are retained.

`MaterialCircuit::new_with_dielectric_capacitors` accepts these material
bindings alongside other authored charge curves, original material resistors,
contacts, mutual inductance and source ports. The selected graph placeholders
are lowered from the material's zero-field differential capacitance. Retained
material/geometry records share immutable source and curve data through native
resistor, contact, conductivity, thermal and magnetic rebuilds. An authored
curve cannot independently claim, alter or remove a source-owned branch.

`rebind_dielectric_capacitors_candidate` explicitly replaces the complete source
set at unchanged geometry and original charge/flux coordinates. Other authored
curves, material records and mutual fields remain intact. The original native
curve-rebind account integrates `V_new-V_old` on the union of both knot sets;
its signed result is parameter work supplied to or returned from the physical
parameter owner, never resistor heat. Ideal-source inconsistency, unsupported
charge, changed geometry, source conflicts and cancellation publish no state.

`tests/material_dielectric_capacitance.rs` covers independent geometry/energy
and bipolar voltage-charge oracles, discharge across an original nonlinear
knot with changing material resistance, source retention through rebuilds,
bidirectional fixed-charge parameter work, original extra capacitor and mutual
magnetic history, source-authority/clamp refusals and cancellation. These
authored fixtures validate numerical behavior, not measured dielectric data.

## Nonlinear instantaneous dielectric response with native polarization memory

`DebyeCircuit::new_with_reversible_instantaneous` combines an explicitly
instantaneous, source-backed D_inf(E) law with the existing source-derived
Debye polarization arms. The original descriptor, RC coefficients, branch
indices, charge/flux ordering and physical resistor-heat owner are reused;
ordinary steps retain the same immutable preparation. The original linear
constructor remains supported. Both profiles use the existing finite thermal
owner within their admitted constitutive scope.

At fixed area A and gap d, the physical coordinates are
q_inf=A*(D-sum P_k) and q_k=A*P_k. The native branch storage therefore realizes
H=A*d*[w_inf(D-sum P_k)+sum P_k^2/(2*chi_k)], with w_inf=integral E dD_inf.
The existing R_k=tau_k/C_k, C_k=chi_k*A/d arms give
tau_k*Pdot_k+P_k=chi_k*E and physical polarization heat
A*d*sum(tau_k/chi_k)*Pdot_k^2. Total terminal charge includes instantaneous and
all retained polarization charges once; it is not the instantaneous charge.
Neither an equilibrium capacitance nor a frequency-dependent C owns this
stored energy. Other original resistor heat stays separately identifiable.

Both source owners must retain the same original card, fixed geometry and
exact non-field context. The nonlinear law is explicitly assigned the
instantaneous role, and its zero-field tangent must match the original
epsilon_infinity to arithmetic roundoff. A total equilibrium D(E) curve is
not an instantaneous source and would count polarization twice. The retained
instantaneous law is field-only; the original background and polarization
increments must be constant-source properties in this closed profile.
Relaxation times retain their original source and nonthermal context; explicit
temperature refresh follows supported positive source paths. Actual nonlinear
voltage, energy and stable energy increments use the original signed native
curve; the original linear source admittance describes only its declared
low-field baseline, not finite-amplitude nonlinear frequency response.

Every accepted step checks both the complete native circuit balance and the
dielectric's own storage-plus-polarization-heat balance against actual terminal
work. Stable constitutive increments and compensated charge transfer retain
small work below a large stationary stored-energy bias. Native capacitor
efforts close each original store's constitutive energy increment; nonlinear
work is not redistributed into another changing capacitor. The independent
dielectric and complete circuit audits retain their original physical
tolerances, and a failed balance is never called heat. Source knots and
relaxation times still need temporal resolution and refinement. Original
source support, native consistent initialization and endpoint constraints,
bounded work and cancellation are admission gates, with no charge projection
or partial memory publication.

The legacy linear source-rebind path refuses the nonlinear profile. Explicit
reversible-source rebinding retains the complete instantaneous law and every
polarization charge and reports the actual fixed-charge parameter work. The
legacy finite thermal entry point below supports temperature-independent
field and polarization potentials with temperature-dependent relaxation times.
The explicit reversible-storage sibling admits the sourced thermal potential
profile described below; no linear C_infinity substitution supplies its energy
or entropy.
This extension adds no saturation of the relaxing polarization, hysteresis,
fractional dispersion, fringing, breakdown or geometry motion.
`tests/dielectric_nonlinear.rs` covers linear and zero-pole limits, nonlinear
open-circuit absorption across a source knot with independent ODE refinement
and both polarities, analytic held-voltage polarization heat, driven warm
memory beside a large independent bias, independent local work closure beside
another changing store, and source/background/context/geometry/clamp/refresh/
cancellation boundaries.
All constitutive data are explicitly authored numerical fixtures.


## Loaded thermal Voce member with physical rate dependence

ThermoplasticSolid composes selected affine-Q Voce hardening with linear or power-law Perzyna kinetics in either selection order. Each physical interval supplies one actual dt to the shared constitutive return; thermal substeps and fixed-point retries retain the incoming history. The six-tetrahedron held-member regression in tests/solid_plastic_enthalpy/voce.rs independently solves axial flow and caloric balance at n=1 and n=3, checks the full storage/physical/numerical energy split, and obtains the same accepted state with changed thermal substeps and enthalpy relaxation. Heating changes tension and lateral shape through the original spatial Newton and enthalpy owners. This is a bounded small-strain synthetic integration check, not empirical calibration, live-runtime cancellation qualification or finite-deformation plasticity.

## Thermal noncondensable gas in native cavitating ducts

`CavitatingDuctPockets::new_with_caloric_gas` binds original
`ClosedGasParcel` owners to the existing bounded DG pocket sites.
`advance_with_gas_heat` advances their mass-preserving internal-energy
increments and supplied heat alongside the SAME duct, cavity volumes,
pressure-work coordinates and collapse events. At each stage original
NASA-9 calorics give T, Cv and p_g=mRT/V; gas storage obeys
dU/dt=heat_power-p_g Q. The shared pressure-rate solve uses
gamma(T)p_g/V and R*heat_power/(Cv V), so temperature-dependent stiffness
and heat-driven liquid motion remain on one trajectory.

Prepared caloric branches amortize continuous-source admission and avoid
constructing full transport receipts at every trial stage. Accepted endpoints
are reconstructed through the original complete material owner. Independent
gas work is obtained from native DG port work plus the reference-minus-vapor
pressure volume term. Both each gas and the combined duct/pocket system must
close their energy balance, using directly evaluated caloric increments rather
than subtracting large absolute formation-energy references.

Gas masses, compositions and source axes remain owned by their original
parcels. Existing vapor and polytropic sites retain their original paths,
including mixed-site collapse events. Domain exits, source refusal, failed
energy gates and final cancellation discard all trial gas, liquid and clock
updates. The heat callback returns signed external watts and must be
deterministic because rejected stages may revisit the same state.

Nineteen native duct/pocket regressions pass: fifteen existing cases plus four
new tests for adiabatic/polytropic equivalence, an independent exact isobaric
moving-column solution with variable Cv and a large formation datum, mixed
caloric-gas/vapor collapse, and atomic late cancellation or final source refusal.
The prepared caloric sampler separately passes all three native tests.
Liquid temperature and saturation pressure remain frozen in this slice.
Vapor phase mass, latent heat, gas dissolution and gas transport must be
supplied by their corresponding physical owners; cavity volume does not
invent those inventories.

## Conservative nonlinear dielectric source replacement

`DebyeCircuit::rebind_reversible_sources_candidate` and its complete-resistor
variant replace an explicitly supplied instantaneous D(E) source and original
Debye source together at fixed geometry, native branch identities and ordered
polarization roles. The proposed pair passes the same source-card/context,
zero-field tangent and temperature-independent-storage admission as initial
construction. Every instantaneous/polarization charge and original flux linkage
is retained bit-for-bit. Source/transformer topology, unrelated storage and
complete original conductor records remain owned by the same native circuit.

The existing capacitor work owner integrates `V_new - V_old` on the union of
the exact original charge knots; each polarization contribution uses its
fixed-charge capacitance-change identity. Compensated summation retains signed
work beneath unrelated stored-energy bias. Relaxation-time changes contribute
exactly zero stored work. This reversible parameter work belongs to the
physical source-change owner and is never deposited as irreversible heat.
Incompatible clamps, insufficient support, switched pole order, missing source
partners, incomplete conductor bindings and cancellation refuse atomically.

`MaterialDielectricCapacitor::at_temperature_query_point` re-queries the same
pinned D(E) claim and complete original field support at a temperature-only
applicability change, retaining every nonthermal coordinate. Its native curve
is shared unchanged, so kinetic tau(T) and conductor refresh can feed the
thermal owner without inventing nonlinear field entropy or storage work. This
does not admit temperature-dependent nonlinear storage, moving geometry,
source-switching impulses or relaxation-model identification. Three focused
tests in `dielectric_nonlinear_rebind.rs` cover independent nonlinear/pole work,
both polarities, large unrelated bias, complete source/history retention,
zero-work kinetic refresh, native-curve sharing and atomic refusals.

## Finite thermal feedback for nonlinear dielectric memory

`circuit_enthalpy::dielectric::step_candidate` admits the nonlinear
instantaneous profile through the existing `DebyeThermalBinding`, finite
caloric network and Newton residual. The original D_inf(E) potential and every
polarization storage increment must be temperature-independent; positive
source-resolved tau(T) and original conductor R(T) may vary. Each trial refreshes
both original dielectric source owners at the same proposed temperature while
retaining every accepted instantaneous/polarization charge, flux, geometry,
source claim and nonthermal coordinate. The unchanged native nonlinear curve
is shared, and the explicit reversible-source rebind must report exactly zero
parameter work for this thermal profile.

The circuit's actual nonlinear instantaneous-plus-polarization energy supplies
field psi=U. Only after the storage-independence gates pass are its field entropy
and fixed-charge field heat-capacity contribution set to zero. The declared
linear low-field background is never substituted for nonlinear energy.
The thermal-offset and reversible-heat contributions are exactly zero; actual
native polarization loss enters the original dielectric caloric cell once,
and each original conductor's loss enters its own receiver. The existing
enclosing account checks the stable electrical internal-energy change plus
finite caloric storage against electrical supply and actual external heat.

The original nonlinear field claim must support the complete declared thermal
interval and charge envelope, in addition to the original caloric/rate sources.
A complete uniform fully solid dielectric volume, fixed density/geometry and
positive combined capacity remain required. This entry point retains its
temperature-independent storage scope. The explicit reversible-storage sibling
below supplies the admitted temperature-dependent field and polarization
potentials. Field support, source identity, ideal-clamp constraints, thermal
event rules, bounded work and cancellation remain atomic gates on both histories.
Endpoint material coefficients make the general thermal feedback first order;
the local energy identity alone supplies no temporal-error bound.

Five groups in `tests/dielectric_nonlinear_thermal.rs` exercise source-knot
crossing and finite self-heating against an independent discrete solve for both
charge polarities; tau(T) history changes and refinement against independent
continuous memory dynamics; external heating with zero poles and unchanged
nonlinear field energy; distinct finite dielectric/conductor receivers with
R(T); and original temperature-support, charge-envelope and late-cancellation
refusals. The accepted local-work case in `tests/dielectric_nonlinear.rs`
separately checks dielectric and unrelated capacitor energy and heat at the
original tight tolerances. All source data are authored numerical fixtures.

## Temperature-dependent nonlinear dielectric energy and finite heat

`ThermalDebyeStorage` retains the complete original nonlinear thermal field
law together with the existing affine Debye storage/rate law.
`DebyeCircuit::new_with_thermal_reversible_instantaneous` admits the explicitly
sourced D_inf(E,T)=s(T)D_ref(E) profile, with positive affine s on a complete
temperature interval. The source scale and reference curve remain distinct
pinned material claims. The declared epsilon_infinity and its temperature
derivative must match the scaled nonlinear zero-field tangent over that whole
interval. The original affine polarization increments may vary with T; their
potentials and derivatives remain owned by `ThermalDebyeCapacitor`.

The same instantaneous and polarization charges enter the original native
descriptor. Current physical charge knots are q=A*s(T)*D_ref(E), with V=g*E;
reference receipts are never treated as current D values. Typed thermal source
rebinding retains the original authority, geometry, branch identities, all
charges/fluxes and conductor records. The immutable authority is shared through
trials. Existing native union-of-knots work supplies the actual fixed-charge
Helmholtz change, including each original pole's parameter work. Legacy closed
constructors and generic source-rebind paths refuse this typed profile rather
than dropping its thermodynamic authority.

`circuit_enthalpy::dielectric::step_candidate_with_reversible_storage` uses
the same finite caloric network, native electrical step and coupled residual
as the existing entry point; `DebyeThermalBinding` is unchanged. It requires
the circuit's exact retained field/pole laws. With z=D_inf/s, the instantaneous
field contributes psi=A*g*s*w_ref(z), S=A*g*s'*coenergy_ref(z), and U=psi+T*S.
The original thermal Debye law contributes only its polarization coordinates,
so its low-field quadratic instantaneous proxy is never counted. The negative
fixed-charge field capacity is admitted together with the original pole and
finite caloric capacities over the complete charge/temperature envelopes.

Stable signed parameter work and the original material delta(T*S) supply
Q_rev=-W_parameter-delta(T*S). This reversible transfer is distinct from native
polarization and conductor heat, each deposited once in its original receiver.
The enclosing balance checks electrical internal-energy change plus finite
caloric storage against actual electrical supply and external heat. Both
charge polarities use the same even potentials. Nonlinear source-knot crossings
and sign changes of the thermal scale slope preserve those source roles.

The profile requires one complete uniform fully solid dielectric volume,
constant density/geometry, positive combined capacity, full source support and
a strictly positive affine scale. General two-coordinate D(E,T) surfaces,
thermal hysteresis, phase changes, breakdown, changing geometry and ideal-source
charging impulses need their own physical owners. Original support, source,
bounded-work and cancellation checks publish neither history on refusal.
Endpoint-temperature coefficients retain the existing generally first-order
coupling; energy closure alone is not a temporal or entropy-production proof.

The `dielectric_nonlinear_thermal` child targets independently check fixed-charge
heating and moving physical charge knots, both thermal/charge signs, negative
combined-capacity refusal, lossless source-driven reversible cooling/heating,
retained typed source authority, incompatible background derivatives and legacy
path refusal. A separate continuous memory oracle solves caloric-plus-field
internal energy for T and integrates the original polarization dynamics with
s(T), chi(T) and tau(T), then compares native timestep refinement. Existing
temperature-independent nonlinear and linear Debye profiles remain regression
cases. All material data are authored numerical fixtures.


## Complete modal loads and bidirectional reflection

`mm_line::continuous::ContinuousModalLine` now drives and observes every mode at both faces when its source has `ModalEnd::Ports`. Each face retains its own original rho*c/S wave reference. The prepared exact owner independently integrates all incident/reflected energies and retains its current-window coupled mean map.

`new_loaded` accepts a quiescent `ModalAdmittanceLoad` containing the complete original pressure-input/volume-flow-output port-Hamiltonian source and an optional original symmetric PSD direct conductance. Every mouth mode is bound in its declared original radial order. Geometry/radius and static-pressure mismatches refuse; interior and load gas density or composition may differ. The source's units, normalization and physical validity remain its producer's responsibility. This declaration never promotes a fitted matrix to a physical energy law.

The direct matrix enters the shared `with_port_conductance` owner as G D G^T; the dynamic source enters the existing skew interconnection before ONE complete exact exponential. All load coordinates, mutual terms, material-memory storage and original resistance survive. Physical load loss remains separate from reference-wave exchange. A purely conductive source introduces no fictitious state: original mouth forcing slots stay present at zero input while only inlet slots are selected as external wave channels. Warm steps reuse the held-input buffer and all exact-owner scratch; refusal preserves accepted physical history and output accounts.

Focused regressions cover bidirectional coherent pulses with unequal face references, a dynamically loaded cylinder against an independent analytic impedance on 12/24/48 axial elements, live mutual-mode reflection with direct and dynamic losses, exact held-wave subdivision, a matched conductive mouth, and source/capacity/transactional refusal. The source and runtime targets together pass 13 focused native checks (5 new and 8 existing), using actual original material/modal/math/algebra/PH source. This is component validation, not a Cargo/DSR gate or full-band BEM radiation qualification. The load may be an admitted complete physical system; no nonlinear forcing, azimuthal modes, sourced higher-mode losses, arbitrary direct nonsymmetric law or closed-loop frequency-error certificate is supplied.

## Coupled capacitance retained through material circuits

`MaterialCircuit::new_with_capacitance_blocks` accepts the native prepared,
reciprocal capacitance blocks on original capacitor branch coordinates, together
with disjoint scalar charge laws, mutual inductance, material resistors and
series electrical contacts. Its `capacitance_blocks()` accessor exposes the
retained immutable authority. The complete matrix, including either sign of
cross-capacitance, owns its branches; their positive graph C coefficients are
topology placeholders. Preparation is shared through circuit copies and rebuilds.

Ordinary resistance, conductivity, contact and finite-body thermal trials retain
the same blocks and physical charge history. Independent electric/magnetic
storage accounting evaluates the complete native quadratic gradient against
charge/flux increments, including the represented-midpoint correction and all
cross terms. Disjoint nonlinear capacitor curves retain their exact integral
increments. This avoids loss of small capacitor work beside an unchanged large
storage bias and prevents diagonal-only accounting from accepting wrong heat.

Scalar capacitor-law replacement retains unrelated blocks and refuses overlapping
ownership. The existing scalar material-storage rebind refuses block-bearing
circuits: changing a matrix requires an explicit full-matrix parameter-work owner.
A moving scalar capacitor may coexist with an unchanged block, but moving a
block-owned coordinate is refused without a complete geometry-dependent storage
potential. Fixed blocks infer no dielectric temperature dependence, entropy or
mechanical work. Their coefficients remain held throughout the existing thermal
fixed-point trials.

The focused material block tests exercise picofarad common/differential RC modes
beside a held 1 F / 10^8 C bias, contact replacement with temperature-dependent
bulk resistance, disjoint charge-law replacement with explicit parameter work,
and two-body finite thermal feedback. They check physical charge, separate bulk
and contact heat, complete storage work, each thermal energy balance, shared
matrix preparation, ownership refusal and cancellation replay. Existing nonlinear
capacitance, source-backed dielectric, finite enthalpy and dielectric-motion
regressions remain green on the selected native source closure. This gate does
not claim a full-workspace Cargo/DSR run.

## Fixed-charge replacement of complete capacitance blocks

`MaterialCircuit::rebind_capacitance_blocks_candidate` replaces authored
reciprocal C blocks at the same original charge and flux coordinates. It
returns the existing `MaterialStorageCandidate`, including the new circuit,
signed per-branch storage contributions and total reversible parameter work.
The original contacts, material resistors, scalar capacitor laws, dielectric
source records, mutual inductance and magnetic source metadata remain retained.
The same native descriptor rechecks the proposed potential against the original
constraints; an incompatible ideal voltage clamp refuses without an impulse.

Every block retains its original set of capacitor branches. Block-list order
and coordinate order within a block may change. Coefficient comparison in
original branch coordinates detects a physically unchanged C and retains that
block's original prepared inverse and ordering, even when another block changes.
An equivalent permutation therefore has exactly zero parameter work and shares
the original preparation. Adding/removing coordinates or regrouping blocks
requires a separate topology and storage-authority owner and refuses here.
Scalar material, nonlinear curve and source-derived scalar dielectric ownership
cannot overlap a block. Empty-to-empty remains a zero-work update.

For a changed block the admitted represented potential determines the work:
W = q^T (C_after^-1 - C_before^-1) q / 2. Inverse coefficients are subtracted
before charge products and compensated row/total accumulation. The method does
not subtract large endpoint energies or voltages. The same graph-branch mapping
is used for both matrix coordinate orders. Per-branch q_i v_i/2 contributions
can be negative; only their complete sum is coupled field storage, and cross
energy is never added again. Displayed endpoint contributions can coincide
after rounding while their separately evaluated parameter work remains nonzero.

Positive work must be supplied by the physical owner changing the matrix and
negative work returned to it; the operation deposits no Joule heat. It is an
explicit zero-time parameter transaction, not continuous material evolution,
geometry/force coupling, a thermal entropy model, a switching impulse solver,
or a new source fit. Existing graph/coordinate capacities bound allocations
and matrix work. Cancellation and native arithmetic/domain refusal publish no
partial charge, flux, material or contact update.

## Coherent finite-film absorption on native layered thermal cells

The opt-in `layered-wave-heating` feature exposes `LayeredWaveHeating` and
`PreparedLayeredWaveHeating`. They bind the existing `PlanarDielectricLayers`
Maxwell owner to one original `BoundLayeredThermalBody`; no second optical or
thermal solver is introduced. One source-backed `DielectricResponse` per
original material layer must retain its conductivity/capacity material-card
identity and agree on shared fixed query coordinates and quantity conventions.
The original optical state and selected constitutive laws remain accessible.
Matching cards do not establish an experimentally matched specimen.

The original layer constructor divides each optical layer at the original
thermal cell faces with the SAME source response in each subdivision. These
equal-medium interfaces are identities. `try_new_cellwise` instead admits one
source response per original thermal cell, preserving each layer's card and
fixed query context while resolving spatially differing optical states. The
original positive even/odd slab loss modes
integrate absorption in each whole control volume, including the standing-wave
cross terms. No depth quadrature, uniform partition of whole-layer heat,
diagonal source-power approximation, or incoming-minus-outgoing heat estimate
is used. Ordered cell IDs, axial thickness and physical area come only from the
thermal body. Preparation refuses more than the native 64 optical films before
allocating; it does not silently coarsen a finer thermal grid.

Both face sources use the same angular frequency, transverse wavevector,
polarization and exp(+i omega t) clock. Inputs are peak tangential E [V/m].
Coherent vectors, explicit second-moment factors C=B B*, and explicitly
independent face sources reuse one prepared field/loss operator. Factors retain
their own mode and layer-mode work limits. Cell power queries allocate nothing
and publish their caller-provided output only after all arithmetic and
cancellation checks. These are spectral-line/time-ensemble quantities, not
PSD/Hz or an inferred model of unknown coherence.

`advance_to` deposits area times each normal absorbed power density times the
physical interval exactly once through `LayeredThermalStep::with_cell_energy_j`.
Maxwell already includes the peak-phasor factor 1/2; coupling does not apply it
again. The native implicit thermal solver retains its nonlinear conductivity,
capacity, contacts, original sources/reservoirs, events, material-law admission
and energy acceptance. The endpoint keeps the authored sources unchanged.
Reflected/transmitted energy is reported separately and never becomes heat.
The returned optical account records its additional cell energy; the native
thermal audit also includes the original external heat channels. Any failed
source query, nonlinear trial, wrong timeline binding or cancellation leaves
the accepted thermal history unchanged.

Optical coefficients and illumination statistics are held over each explicitly
requested interval. This is a fixed-geometry, frozen-optical-coefficient
approximation: the base owner retains its supplied optical state. The
`ThermoOpticalWaveHeating` wrapper below provides accepted-temperature source
refresh. Implicit thermo-optic feedback, phase change, thermal emission,
broadband emissivity, mechanical feedback and measured-material accuracy are
not claimed.
Five focused test groups in `tests/layered_wave_heating.rs` cover split/unsplit
film agreement, full source coherence, independent analytical weak-loss cell
integrals, actual native enthalpy/temperature with authored heating, reversal
and area scaling, opaque two-face heating, and transactional source/work/timeline
refusals. Numerical execution is a separate validation step.

### Thermo-optic feedback from the accepted spatial temperature

`ThermoOpticalWaveHeating` binds one original optical state and its immutable
`MaterialCard` per thermal layer. Its sealed generic source defaults to
`VisibleConductorStatePoint`, preserving existing conductor callers. Choosing
`VisibleOpticalStatePoint` admits mixed metal, Cauchy, tabulated glass and
Sellmeier layers through the same solver. Dielectrics retain their original
spectral intensity absorption and phase-index law; RGB-only absorption refuses.
The complete source includes the same-state magnetic permeability; no DC loss,
substitute spectrum or nonmagnetic default is added. Its temperature coordinate and shared fixed
conditions must agree with the original conductivity/capacity owners. Each
step obtains actual cell temperatures from the accepted native enthalpy state,
then refreshes every cell's original selected optical/mu claims at that temperature.
The other query coordinates and source pins remain fixed. No layer-average
temperature or separate optical temperature input replaces the spatial state.

The original coherent film owner treats these cells as distinct physical
sublayers, and the existing implicit thermal owner advances their absorbed
energy once. Before returning, every optical source is readmitted at its
proposed thermal endpoint. If the endpoint leaves the original optical source
support, introduces an unsupported dispersion pole, fails passivity or
representability, or is cancelled, no joint endpoint
is returned and the accepted thermal history is unchanged. Results expose the
applied and endpoint source snapshots in original cell order; refreshing the
endpoint is not a second heat credit. The next interval uses its accepted
temperature and therefore closes the explicit optical-to-heat-to-optical loop.

Optics are held during each requested interval. This is a fixed-geometry,
quasi-steady explicit split; algebraic heat acceptance does not certify its
temporal accuracy. Callers retain step-refinement responsibility. The wrapper
does not infer implicit optical tangents, transient electromagnetic storage,
geometry-change work, broadband thermal emissivity or a phase-change law.
Focused tests compare nonuniform feedback with direct source refresh and an
independent two-cell caloric solution, exercise successive accepted steps,
and verify late-cancellation, source-substitution and endpoint-domain rollback.
Mixed glass/metal tests retain each spatial source family and local absorption
through successive native heat steps, including atomic refusal of a glass
endpoint outside its original optical source support.

## Coherent porous acoustic absorption on native thermal cells

`porous_heating::PorousHeating` binds one original
`fs_duct::porous::PorousLayerExcitation` to a complete ordered depth partition
of existing `fs_thermal::enthalpy` cells. Each distinct CellId receives its
analytic viscous and thermal-relaxation work over its own interval, including
coherent standing-wave cross terms. Heat is not inferred by subtracting face
powers or divided uniformly over the specimen. The original peak-phasor 1/2
and full cross-sectional area are already included and are never applied again.

The mapping admits adjoining, straight, uniform prismatic control volumes:
cell centreline endpoints are opposite face centres, increasing depth follows
the original acoustic layer normal, volume equals full face area times depth,
and the complete partition ends at the original layer thickness. Length,
alignment and volume comparisons allow floating-point roundoff only. Duplicate
destinations, incomplete coverage, gaps, bends, wrong area/volume, inactive cells
and stale geometry/mass/material bindings refuse. Material storage basis and
available original caloric identity remain bound. A caller-declared preparation
budget is capped at 4096 cells before allocation. One analytic integral per
cell prepares the heat; thermal substeps only copy the cached source map.

`advance_to` uses the existing transactional `DepositionObserver::sources`
path. The original explicit thermal timeline retains its conduction CFL,
native caloric domains and latent storage, authored heat, contact/reservoir
exchange, events and acceptance gates. Generated acoustic watts are added for
the current substep and never persisted in authored boundaries. The caller
must not also author the same acoustic heat as a boundary source. Exact native
step durations accumulate per-cell heat, both dissipation mechanisms and the
incident/outgoing acoustic energy in compensated accounts. Only cell work is
deposited; acoustic boundary closure is checked independently without scaling
the heat. The native thermal report already includes that heat together with
its other external exchanges. Failure or cancellation publishes no candidate.

The caller explicitly authors each cell as one thermally lumped specimen
volume and owns the justification for one fluid/skeleton storage temperature.
The viscous/thermal acoustic mechanism split does not determine constituent
heat fractions. Acoustic coefficients and a monochromatic coherent drive are
held fixed over an advance: temperature-dependent acoustic refresh, frame
motion, startup transients and two-temperature heat exchange are not inferred.
This composes the original rigid-frame constitutive and thermal models and
does not extend their source support or provide material qualification.

Five focused regressions in `tests/porous_heating.rs` cover independent
standing-wave volume integration and 1/3/11-cell partition invariance, native
finite storage, original conduction/reservoir/authored heat, scheduled boundary
changes and continuation, latent storage and domain/cancellation rollback,
exact native zero-drive behavior and physical-mapping refusals.

## Adhesive snap release on the original finite caloric cells

`friction_heat::adhesive::snap::JkrSnapHeat` retains a sealed original JKR
branch and regenerates every accepted monotone path segment before assigning
its jump release. The expected start version and source-bound final state must
match; a stale retry or substituted source cannot silently deposit the same
accepted release again. The trajectory owner must supply every path reversal.
Smooth reversible normal work deposits exactly zero. Exact attachment and displacement-spinodal releases
remain with `SourceJkrSphereContact` rather than caller-editable jump reports.

An explicit `fs-tribo::HeatPartition` assigns surface A, surface B and named
outside work. `DissipationStep::from_work` accepts the original integrated
joules without inventing a power or duration; the existing `WorkLedger` retains
their channel totals. The two ordered normal material-card identities must
match the original `EquilibriumEnthalpyPhaseCurve` objects actually registered
in the native timeline. Cell, segment, body, material, mass and geometry remain
pinned. Heat-only mapping adds no conducting link or fictitious solid-contact
area, including at detachment.

One original `DepositionSimulation::apply_heat_j` updates the finite caloric
state for a batch. Its per-cell and summed absolute storage gates, a separate
release-transfer gate and the final cancellation check admit the thermal and
branch-history candidates together. Invalid paths, changed receivers, material
domain failure, exhausted budgets or cancellation publish neither owner.
The native thermal audit labels this impulse external to the thermal subsystem;
it is an internal debit/credit in the contact/thermal/outside work account.

This is a same-time, frozen-normal-source heat split. A trajectory owner
localizes jumps and supplies every monotone segment; thermal transport between
events requires separate native timeline advances. The adapter does not infer
heat fractions, solve radiated waves, refresh temperature-dependent normal
coefficients, or claim whole-system internal-energy closure from enthalpy alone.
Seven focused groups in `tests/adhesive_snap_heat.rs` exercise the original latent
chart and heat update, a complete adhesive cycle including exact-fold release,
zero heat on smooth tensile/compressive paths, source/history/material refusal,
every cancellation phase with deterministic retry, strict transfer closure,
and integrated-work units.

## Finite-time viscous-support heat on the original thermal timeline

`JkrSnapHeat::propose_timed_damped_normal_motion` consumes the original
`JkrNormalMotionTrace` and a separate `JkrViscousHeatRoute`. The route
explicitly selects surface shares for the two already-bound caloric cells and
names any outside destination. It never inherits the adhesive snap partition.
Sealed support work and smooth-segment timing are distinct from original snap
release, reversible contact work and spring storage.

A private observer adds each original segment's average dissipative power
through native `DepositionObserver::step_sources`. The original thermal
solver clips at power changes and transports/cools heat while it is produced;
authored sources, boundary schedules, reservoirs and thermal topology retain
their original owners. All-outside and zero-loss routes add no thermal source
or segment clipping. Original snap events still split the same timeline and
retain the native scheduled-event ordering at each timestamp.

The accepted native thermal durations determine actual injected energy.
Per-segment, aggregate and cumulative absolute transfer gates compare that
energy with original mechanical work; clock rounding cannot be repaired by a
manufactured endpoint heat impulse. Relative native transport/snap audits and
support-transfer residuals compose with compensated accumulation. Reports
keep both work ledgers and outside destinations separate and identify authored
transport heat after removing the internal support source. Full continuation,
receiver identity, the shared transport substep allowance and atomic refusal
remain binding.

The reconstruction uses the original mechanical midpoint loss and constant
power within each accepted mechanical segment. It therefore needs mechanical
and thermal timestep refinement; it is not a resolved within-step viscous
waveform or a temperature-dependent contact law. Normal source coefficients
and both explicitly declared heat partitions remain frozen. The implementation
reuses the existing thermal solver and common event-timed snap engine.

## Event-timed normal snap heat through original transport

`JkrSnapHeat::propose_timed_normal_motion` advances the original thermal
timeline to every sealed normal-motion event, deposits its regenerated release,
then continues transport to the interval endpoint. The explicit clock origin
maps mechanical time into thermal time; the initial clocks and full motion
continuation must match. Nonzero event intervals must remain representable on
that clock. Native scheduled births, boundary changes and thermal contact
events at a snap timestamp execute before its heat impulse, in their original
order. Normal snaps do not invent conducting area or thermal links.

Transport uses the original `advance_observed_limited` no-op-observer path;
each release uses the original `apply_heat_j`. The shared substep allowance is
passed as the remaining native work limit, so splitting at events does not
multiply the caller's transport budget. Original finite calorics, receiving
cell/material identity, geometric/mass binding and heat partition remain in
force. Complete candidate heat history and thermal state publish only after
all interval, transfer, cumulative accuracy and cancellation gates succeed.

Reports retain each physical timestamp and native transport/impulse audit,
separate ordinary boundary/source heat and incoming material energy from the
internal snap transfer, and accumulate native relative-storage residuals with
double-double arithmetic. Large absolute caloric reference energies do not
replace those small-work gates. Full motion chronology still prevents replay
with a fresh heat stamp after a closed contact cycle.

The normal source and heat fractions remain frozen; event-timed one-way heat
coupling does not imply temperature-dependent adhesion or resolve contact
waves. The existing same-time proposal remains available with its own declared
scope. Four focused groups cover original split-timeline equivalence and a
material difference from endpoint-lumped heat, exact scheduled-event ordering,
clock mapping, ordinary zero-snap transport, shared work limits and atomic
caloric/transport/cancellation refusal.

## Adhesive body islands with source-resolved finite thermal transport

`friction_heat::adhesive::propose_rigid_adhesive_island_heat` couples original
JKR/DMT/Maugis solid-contact disks to the existing simultaneous native body
island and finite enthalpy network. Each attached patch has its own solved face
temperature, original absolute-temperature axis and selected junction-strength
claims. Zero or tensile net normal load does not remove its real solid area.
Detached disks cannot acquire a conducting face through this entry; the normal
and topology owners must first remove or replace that physical thermal route.

The prepared full body response is reused at every thermal trial. Canonical
mechanical IDs and independent thermal-link IDs preserve the mapping into shared
bodies and cells. Actual combined sliding/spin dissipation enters one heat owner
per patch. The original series-conduction link is replaced within the private
thermal trial, so conduction is counted once. Per-side endpoint conductivity
and original geometry determine k*A/d; the thermal face must match the actual
solid disk's centre, normal and area. Signed heat may leave a hot body even while
the adhesive junction dissipates. No prescribed heat fraction or added thermal
store supplies the interface temperature.

Complete native impulses, source receipts and caloric endpoints are returned
only after local work gates, the original thermal gates, cancellation, and both
an independent `delta K + delta U - outside heat` check and summed audit check.
Invalid source support, route identity/geometry, work budgets or late cancellation
leave accepted state unchanged. The heat callback owns other transport only.

This is explicitly a frozen-normal and frozen-geometry thermal split: junction
strengths refresh, while normal elastic/surface-energy coefficients, contact area
and poses remain with their existing owners. The perfect-face orthogonal-cell
backward-Euler approximation needs spatial/time refinement; convergence does not
establish uniqueness or validation. Event-resolved stopping, thermal expansion,
normal thermal feedback and evolving attachment topology are not supplied here.
Four G1/G3/G4 groups in `tests/adhesive_friction_heat.rs` compare shared-inertia
pure-spin work and a two-face/three-cell heat solve with independent analytical
elimination, including zero/tensile load, source feedback, signed ordinary
conduction, permutation and domain/geometry/cancellation rollback. Prismatic
face-area admission tolerates only arithmetic roundoff in reconstructed volume
over length; a separate decimal-coordinate case still refuses physical excess.

### Adhesive trajectories through arrest with finite thermal feedback

`friction_heat::adhesive::transient` closes one original rigid adhesive pair's
complete sliding/spinning trajectory against the existing finite enthalpy
network. Geometry, complete native mobility and normal-response rows are
prepared once. Each thermal trial refreshes the same pinned junction claims at
its solved face temperature and advances the original trajectory from the same
accepted momenta. Curved motion, the exact-ray fast path and finite extinction
retain their original work, stop-interval and unresolved-travel semantics.

Only original integrated traction power and the physical extinction-tail loss
enter the massless-face balance. Each bulk retains its own finite caloric state
and endpoint conductivity; signed cooling and ordinary conduction continue after
mechanical arrest. The original contact link is replaced once within the private
trial, without adding a heat store, prescribed heat fraction or duplicated loss.
Actual kinetic-plus-caloric change and the summed independent audits must both
close. Preparation is immutable; the apply entry publishes both body momenta,
the thermal network and staged source payload only after final cancellation.

Strengths are sampled at the thermal endpoint and held over the mechanical
phase. Temporal subdivision is required for material feedback and heat-release
timing; the mechanical trajectory's error estimates do not bound this thermal
sampling error. Normal coefficients, solid area, poses and inertia remain fixed.
Normal thermal potential, expansion, moving footprints and attachment topology
are not supplied by this phase. Four G1/G3/G4 groups in
`tests/adhesive_friction_heat/transient.rs` compare finite-inertia spin arrest and
thermal strength feedback with independently eliminated cell balances, preserve
the original curved native trajectory, and exercise post-arrest conduction plus
source/budget/late-cancellation rollback and deterministic retry.

## Source-bound wear on the actual directional friction and heat trajectory

`friction_heat::directional::wear::apply_source_directional_trajectory_heat_with_wear`
advances original directional rigid recoil, source history, finite enthalpy and
both surface Archard volumes in one transaction. The existing full-matrix
curved sliding/arrest owner supplies actual unsigned travel, not net displacement
or an impulse-derived distance. The accepted original friction query supplies
wear coefficients, and each original bulk hardness is refreshed at its own
actual caloric temperature.

Body and thermal-cell budgets are checked before staging. Every source/wear
check and final cancellation precedes joint publication; a late hardness-domain
refusal leaves bodies, history, thermal state and both volumes unchanged. Wear
creates no energy share or replacement heat partition, preserving the original
mechanical/thermal energy account. Geometry, poses and normal load retain the
original phase assumptions. Endpoint coefficient/hardness sampling requires
time refinement; predicted volume does not infer geometry recession, debris or
mass removal. Native roundoff slip is retained rather than silently clamped.
Four physical test groups cover curved travel and distinct thermal hardness,
forward/reverse sliding and arrest, zero wear, original source/role/fixed-context
admission, unchanged heat and late atomic rollback.

### Heterogeneous porous stacks on the same thermal timeline

`porous_heating::PorousStackHeating` binds one complete local-depth thermal
partition to every layer of an original `PorousStackExcitation`. The original
scattering solution supplies each layer's coherent internal drive, so interface
reflection, transmission and two-sided interference determine the heat in each
material. Applying the external drive independently to each slab is not used.
Only the stack's external incident/outgoing waves enter its acoustic energy
account; power exchanged between adjacent layers is internal, not extra heat.

Preparation reuses the single-layer analytic cell integrals and physical
admission, then merges all sources into one cached native heat observer. The
complete specimen shares the 4096-cell ceiling. Cross-layer cells must be
distinct and adjoin on one straight, full-area body normal. Each material's
original mass, caloric/storage law and identity remain bound. Authored thermal
contacts across interfaces participate in the same native advance, including
different heat capacities or latent storage. Cell depths stay local to their
authored acoustic layer; `layer_cells` exposes that exact partition.

The existing fixed-drive and thermally lumped fluid/skeleton assumptions remain
unchanged. Four focused stack regressions cover heterogeneous internal-field
volume integration and partition invariance, distinct native material storage,
cross-material conduction with authored heat and events, one-layer bit parity,
continuation/cancellation, and complete physical mapping. They establish the
adapter's numerical behavior without extending material source validity.

## Driven Hertz sliding, spin and rolling with finite thermal feedback

`friction_heat::hertz_rolling::propose_rigid_hertz_rolling_driven_heat` admits
held world forces and body-frame COM torques through the original five-coordinate
native contact response. The coupled thermal root resolves the selected surface
friction at actual interface temperature and admits bulk rolling at the original
sphere temperature. Surface heat enters the original common face; rolling heat
enters the original whole-sphere caloric cell. Full-body actuator work remains a
separate mechanical supply in `delta_K + delta_H - Q_external - W_drive`.

The corresponding `collinear::propose_rigid_hertz_rolling_collinear_driven_heat`
uses the existing continuous invariant-subspace trajectory. Static retention,
breakaway, arrest and reversal keep their original event times and opaque
continuation. Each thermal trial starts from accepted bodies/history, and exact
phase-integrated surface/bulk losses retain the same caloric destinations.
Common-body acceleration changes actuator work without creating friction heat.
Both original unforced APIs delegate with zero wrenches.

The full five-coordinate branch retains its frozen-geometry held-effort model
and unresolved-event refusals. The collinear branch retains its admitted
invariant subspace and held coefficients/load. Neither adds pose/normal
integration or changes heat partition. Endpoint thermal sampling still requires
time refinement. Seven new physical groups cover an independent forced
mechanical/thermal reference, mixed shear/spin/rolling, driven reversal and
thermal weakening, common-acceleration/static limits, exact zero-drive behavior
and late atomic source/cancellation refusal; eleven existing tests also pass.

## Powered shared-body Hertz assemblies with finite thermal feedback

`friction_heat::hertz_rolling::multiple::propose_rigid_hertz_rolling_island_driven_heat`
accepts held external world forces and body-frame COM torques in complete
original body order. The prepared native island retains full shared-body
coupling and applies each participating body's drive once. Original per-face
source temperatures and conductivities still determine each interface's surface
heat; every sphere's patch bulk losses accumulate once into its unique original
caloric inventory.

`external_mechanical_work_j` remains an independently audited supply in
`delta_K_unique_bodies + delta_H - Q_external - W_drive`, including common-body
acceleration. All source and thermal trials restart accepted bodies and loaded
rolling histories, with one final immutable proposal after cancellation and
energy checks. The original unforced entry retains its native constructor and
behavior without allocating a new whole-body zero-wrench input.

Three new physical groups compare shared-sphere force/torque/three-capacity
response with independent equations and temperature-dependent sources, check
exact zero drive and common acceleration without additional physical heat, and
exercise late source/cancellation refusal from loaded histories. All 27 shared
Hertz/rolling finite-heat tests pass. Native spectator/normal-drive refusals,
fixed geometry/load, held-effort event limits, original source applicability and
thermal time-refinement requirements remain.

### Native temperature feedback into coherent porous absorption

`porous_heating::feedback::TemperaturePorousHeating` binds each original thermal
cell to an unchanged `PorousFluidLaw` and treats its physical length as one
homogeneous acoustic sublayer. `IsobaricPorousGas` retains the original `GasSpec`,
prescribed pressure and explicit temperature support; the native `GasState`
constructor admits the support endpoints and each evaluated temperature. Native
gas density, sound speed, viscosity, conductivity and pore drag therefore evolve
together. The original thermal cell supplies its actual caloric temperature,
geometry, mass and storage law; temperature-dependent pore density never
replaces the specimen's thermal mass or caloric identity.

Only changed cell temperatures reevaluate their native gas and porous
coefficients. Any changed layer triggers one shared coherent stack solve, so
temperature gradients alter interface reflection and the internal drive of all
layers. Unchanged temperatures, including a latent plateau, reuse the prepared
analytic sources. Physical admission reuses the complete straight, adjoining,
full-area stack binding; the existing 4096-cell ceiling, composition budget and
cancellation gates bound preparation and refresh work.

`PorousHarmonicDrive` fixes the external peak pressure phasors, frequency, face
area and positive real wave reference. Acoustic source energy counts only the
two exterior faces. Each native thermal step deposits its accepted-start
viscous and thermal watts exactly once through the shared heat observer and
energy account. Endpoint refresh prepares the next rate and validates the
returned cache without changing already deposited energy. Original conduction,
authored heat, boundary events, latent storage and continuation remain native;
late source-domain or cancellation failure publishes neither candidate.

The default model assumes rigid-frame, normal-incidence quasisteady harmonic acoustics;
the opt-in conserved-tangential-wavenumber mode is specified below. Both assume
prescribed isobaric background composition and pressure, and a justified common
gas/skeleton temperature. The authored caloric law must represent storage under
that background. Finite gas inventory, reservoir mass transport, evolving mean
pressure, transient acoustic storage and parametric work are not inferred.
Native conduction CFL alone does not control nonlinear absorption feedback;
callers must refine thermal timesteps against their source dynamics.

Five focused physical regressions cover native nonuniform gas-state and
coherent-loss parity, latent/cache behavior, events and continuation, atomic
domain/mapping/cancellation refusal, and accepted-start rates plus timestep
refinement against an independent two-cell direct-owner RK4 reference with
conduction, authored heat and outgoing acoustic energy. The nine existing
fixed-source porous thermal regressions remain unchanged and pass.

### Original liquid-source feedback in the same porous heat owners

The harmonic, spectrum and covariance heating owners accept a defaulted
`PorousTemperatureFluid` source; existing `IsobaricPorousGas` constructors and
gas-state views retain their original behavior. `IsobaricPorousLiquid` supplies
one concrete liquid-saturated path through the existing
`ResolvedLiquidDensityState`. Its original selected `rho(T,p)` surface and Cp
derive the same state's thermal expansion, compressibility, Cv and adiabatic
sound speed. Viscosity and conductivity refresh with them; an original kinematic
viscosity is converted with the new density. No ideal-gas closure, fitted
temperature multiplier or alternate porous response is introduced. The native
local density response retains its own single-phase, derivative-support and
thermodynamic consistency limits. The general-fluid use of Cp/Cv and adiabatic
bulk modulus follows the [COMSOL poroacoustics general-fluid derivation](https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.182.html).

The liquid binding retains the exact immutable card, typed temperature/pressure
axes, selected claims and other original context, including composition. It
admits explicit thermal-support endpoints and every later native temperature
through the original continuation API; endpoint admission does not fill source
holes or establish intermediate phase stability. `LiquidSource` refusals retain
the native liquid error, original card identity, complete reference query and
attempted temperature. Source and cancellation failures publish neither the
candidate acoustic cache nor the candidate thermal state.

This is isobaric constitutive continuation of a thermally lumped rigid specimen.
The original thermal law must justify common liquid/skeleton temperature and
storage under the prescribed background. Density changes do not move fluid,
change pore inventory, add stored mass, or alter fixed specimen geometry and
exterior wave references. Native authored masses and caloric identities remain
unchanged. Actual positive per-layer viscous/thermal integrals still supply the
accepted-start thermal watts; external power residuals never become heat.

One original liquid evaluation serves each changed physical cell across all
frequencies/angular modes. Immutable evidence-bearing states are shared across
spectral caches and candidate clones. Unchanged temperatures and latent plateaus
reuse them. `PorousFeedbackWork::fluid_evaluations` counts all changed fluid
evaluations; its retained `gas_evaluations` field counts only gas sources.
Focused regressions cover direct original-stack parity at oblique incidence,
local density/Cp/Cv identities, native accepted-start energy and unchanged mass,
zero thermal-expansion loss, shared spectral/covariance source work, and atomic
source-context/cancellation refusal.

<!-- porous-spectrum-feedback:start -->
### Stationary finite-band porous heating with native temperature feedback

`porous_heating::feedback::spectrum::PorousSpectrumHeating` advances a common
physical specimen under distinct statistically independent frequency bins.
Each bin uses the original coherent two-sided porous stack; its positive
viscous and thermal cell watts add to one native heat source. Complex pressures
are never summed across frequencies. Bins are canonicalized by increasing
positive angular frequency. Duplicate frequencies refuse: coherent waves at
the same frequency must be combined in one authored bin. All bins share the
same full specimen area and original thermal/material bindings.

Each `PorousHarmonicDrive` contains the bin's already integrated peak pressure
phasors and a fixed positive real exterior reference. The explicit
`angular_pressure_density_bin` helper accepts one-sided angular-frequency
amplitude density, with mean-square pressure `integral S_pp(omega) d omega`,
and a positive quadrature weight. It converts density to peak pressure using
`sqrt(2 * weight)`; bandwidth and the harmonic factor are applied exactly once.
The two ports remain fully coherent within that statistical component.

Changed native cell temperatures evaluate the original gas law once per cell,
shared across all frequencies. Each affected bin rebuilds its changed porous
coefficients and performs one coherent solve. Unchanged/latent temperatures
reuse prepared rates. Positive powers and exterior accounts use compensated
sums; per-bin energy is diagnostic, never an extra thermal credit. One shared
accepted-start observer preserves original conduction, authored heat, events,
latent storage, continuation and atomic endpoint/cancellation refusal.

Explicit bin and cell-times-bin budgets precede allocation, with hard ceilings
of 4096 bins and 65536 cached cell/bin pairs. Existing native composition limits
and the 4096-cell heat ceiling also apply. The model retains the original
isobaric background and thermally lumped storage assumptions and represents
stationary ensemble or long-time means. Finite-time beats, acoustic transients,
frequency conversion and within-bin partial coherence require other owners.
Spectral quadrature and thermal source feedback each require refinement.

Five focused groups cover independent spatial work, within-bin phase effects,
angular-density normalization and spectral quadrature refinement, shared caches
and latent storage, exact single-bin/zero-drive behavior, events/continuation
and atomic refusals, plus two-frequency temperature/energy refinement against
an independent direct-owner RK4 reference. Existing fixed and single-frequency
thermal tests retain their original source and arithmetic.
<!-- porous-spectrum-feedback:end -->

## Prescribed source-frequency sweeps in modal performances

The shared physical-tone schedule accepts optional `frequency_sweeps COUNT`,
`sweep TONE_INDEX POINT_COUNT MAX_FORCE_ERROR_N`, and `omega SAMPLE RAD_S`
records after gain envelopes. Indices refer to original physical tone rows.
Knots begin at sample zero, increase strictly within the declared horizon,
remain below the original guarded Nyquist limit, interpolate linearly, and
hold the final frequency. A zero target frequency is a stopped prescribed
phase, not a change to resonator stiffness or a string-tension control.

The schedule admits a conservative whole-curve bound on complex-envelope
interpolation error scaled by the original force quadrature magnitude in
newtons. For a sample of duration h, fixed carrier wc, gain extrema G and
maximum gain change per sample Dg, detuning D=max|w-wc|h, and maximum frequency
change per sample Dw, the bound is |F|[2 Dg D + G(Dw h + D²)]/8. It is checked
against the explicitly declared per-tone force-error allowance. Independent
gain/frequency extrema may overestimate the bound; the report is the largest
admitted per-tone bound, not a total across unrelated ports or a pressure
bound. This floating-point evaluation is not an interval certificate including
accumulated phase roundoff or an assertion of bandlimited modulation sidebands.

Knot and metadata storage is reserved before allocation, with at most 65,536
aggregate frequency knots. Candidate interval queries allocate nothing and do
not advance the accepted clock. Cursors publish only after physical and causal
observation acceptance. Focused tests cover original-row addressing, exact
integrated phase through silence and zero-frequency turns, held tails,
force-unit and signed-gain error admission, malformed knots, storage and
guarded-clock refusals, and repeated/rejected candidate queries.

### Phase-continuous physical source frequency sweeps

The independent broadband performance path accepts optional `frequency_sweeps`,
`sweep TONE_INDEX POINT_COUNT MAX_FORCE_ERROR_N`, and `omega SAMPLE RAD_S` rows
after gain envelopes. Original physical tone indices survive voice grouping.
The shared parser admits the complete target-frequency and gain curves against
the original guarded clock and explicit absolute physical force error allowance.
The carrier frequency used for prepared propagation remains fixed; these are
source gestures, not changes to structural modes or materials.

One immutable modulation plan per swept tone/sample generates both endpoint
quadratures in the same interval-start carrier basis and advances the unscaled
phase by the integrated target frequency. Signed gains, silence, and re-attacks
retain that phase. The existing affine source/receiver owners propagate the
reconstructed force together with physical work and heat. This is a bounded
force reconstruction of a chirp, not an exact chirp solution or an acoustic-error
certificate. Original mechanics, receiver and resource gates remain active.
Sweep state and additional per-tone/per-mode modulation work consume the
existing source budgets. Curves and phase publish only after accepted steps.

The v1 sidecar exposes `physical_force_frequency_sweeps` and
`maximum_tone_force_reconstruction_error_n`; the example is
`modal_swept_tones.performance`. Focused checks compare against an independent
physical chirp ODE, the analytic integrated phase across silence, signed modal
projection, callback equality, source-free ringdown, explicit admission, and
warm rejected-step retry. Empty sweep sections preserve existing envelope bits.

## Coupled performance source-frequency sweeps

The v2 `CoupledTonePerformance` consumes the same original-row frequency
schedule and physical force-error admission as independent performances.
Sweeps opt into the existing complete affine source/receiver preparation;
the prepared carrier frequencies, structural modes and material laws stay
fixed. Different physical sources at one prepared carrier retain separate
phase histories before endpoint coalescing. The exact carrier operators
propagate their reconstructed force, original work and heat channels, and
causal receiver on the same physical interval.

The driver queries each sweep once per sample and reuses its cached proposal
across signed modal projections. Unscaled phase follows the exact target
frequency integral, including silence and re-attack. It stages gain-scaled
endpoints and the next phase before touching accepted held inputs; source,
observer, gain/frequency cursors and held events publish only after acceptance.
A refused physical step leaves the previous state and retry proposal intact.
Release retains coupled motion and receiver/delay history.

`frequency_sweep_count` and `maximum_sweep_force_error_bound_n` expose
the authored count and largest per-original-tone reconstruction bound [N].
The existing CLI writes those values and the reconstruction scope into its
sidecar. The bound concerns the prescribed force reconstruction, separately
from receiver-fit error, accumulated floating-point phase error and modulation
bandwidth. It does not modify string tension or infer a nonlinear excitation.

`examples/coupled_swept_tones.performance` exercises coherent signed forces,
distinct sweeps sharing a carrier, gain releases/re-attacks and the existing
causal delayed receiver. Focused G3 file tests check the independent integrated
phase, bitwise callback partitions, retained state and positive-heat ringdown
after source release, explicit admission, rejected-step retry, and exact
empty-sweep parity. Both existing CLI children have been exercised on their
swept examples through the unchanged decimation and PCM owners.

## Retained Mindlin microslip on shared bodies with finite contact heat

`friction_heat::mindlin::propose_rigid_mindlin_island_heat` and its driven
counterpart couple original source-bound scalar Mindlin histories to the existing
finite multi-interface thermal network. Actual native masses, full inertia and
attachment arms determine shared-body recoil and independent actuator work.
Each physical microslip/gross-sliding loss enters its massless common-temperature
face once; conduction and partition follow original finite-cell geometry and
endpoint conductivities, without a caller heat fraction or added heat store.

The mechanical island is solved once. Each thermal trial resolves the original
bulk and interface sources at the actual separate cell/face temperatures, then
uses exact source-owned temperature continuation to prove unchanged mechanical
coefficients throughout the interval. The same retained spring history and
physical dissipation therefore remain valid through thermal Newton iterations.
Nonconstant mechanical claims refuse; temperature-dependent caloric/transport
laws retain their original owners and the existing numerical solve.

Admission checks original circular footprint/face geometry, native attachments,
initial and endpoint bulk temperatures, solved interface temperatures, unique
heat-link ownership and consistent native-body/full-card ownership of shared
thermal cells. The source callback owns original caloric/conductivity bindings.
Native kinetic energy, retained shear storage, caloric energy, external heat
and once-per-body actuator work close before bodies, histories, refreshed source
laws and the thermal candidate return together. Optional entropy uses the
original network audit. Cancellation or late source/refinement/budget failure
publishes none of the candidate states.

Four G1/G3/G4 groups compare cell and face temperatures with an independent
three-capacity temperature matrix, retain loaded reversal memory over repeated
thermal updates, check signed hot-to-cold flux and external work, canonical
contacts/common acceleration, original source temperatures and atomic refusal.
The fixed-normal/frozen-geometry phase does not supply changing-footprint heat,
thermal expansion or history transfer for varying elasticity/friction; finite
thermal space/time discretization retains its refinement requirements.

## Complete sliding-film layer heat on the impact clock

`render::plate::impact::material::heating::SlidingFilmLayerHeatRoute`
binds every retained lower-to-upper layer of one selected original sliding
channel to an existing finite caloric reservoir. Uniform channels retain
one implicit layer. Film indices use the separate sliding insertion order;
normal-pressure film indices and liquid feedback bindings cannot alias
these routes. Coverage is complete, duplicate physical-channel routes
refuse, repeated reservoir destinations are permitted, and unselected
channels remain mechanical losses outside the chosen thermal boundary.

`HeatedImpactSystem::new_with_sliding_losses` combines these routes with
the existing contact, normal-pressure film, stratified normal-film, and
material-memory routes. Legacy constructors delegate with no sliding
routes. Original wall ports, immutable layer interfaces, viscosity
profiles, and the native force/tangent owners remain authoritative.

The heat stage evaluates each selected sliding island once at the original
midpoint configuration and discrete-gradient velocity. It consumes the
native per-layer integral of `tau^2 / mu`, including the signed
pressure/Couette cross term; pressure-only heat fractions cannot substitute
for that combined physical work. Only the floating-point closure remainder
is assigned to the current largest layer contributor. Each selected channel
counts once against the original mechanical dissipation allowance,
independently of its layer count. The large bounded layer report lives in a
separate non-inlined helper, so normal-film/contact-only calls retain their
existing stack requirements.

Motion, finite enthalpy stores, thermal transport/network state, and the
original clock publish through the existing atomic acceptance transaction.
Cancellation, unrepresentable heat, invalid routing, or a late material-chart
refusal leaves accepted state unchanged. Direct use of the heat-route
constructor holds sliding constitutive profiles fixed. Temperature feedback
requires the explicit liquid-film binding below; wall heat destinations and
the original thermal clock remain unchanged.

`LiquidFilmImpactSystem::new_with_sliding` adds explicit source-temperature
feedback through the SAME finite-store proposal and accepted mechanical clock.
`LiquidSlidingFilmLayerBinding` covers every original layer of every channel in
each selected sliding island; a uniform passage has one implicit layer. Its
temperature reservoir must be its existing heat destination. Normal-pressure
and sliding-film insertion indices remain separate even when both are zero.
The original material-card identity, pinned source claims, fixed nonthermal
context, fully liquid phase and constant density/mass/volume gates apply to
every row. There is no layer inventory, pressure-derived temperature or implicit
wall heat partition.

At the accepted endpoint, original local temperatures requery viscosity and
the native `ShearedFilm::with_layer_viscosities` recompiles hydraulic mobility,
Couette entrainment, harmonic wall drag and complete layer work together.
Geometry, wall rows, normalized interfaces and the common Barus coefficient
remain fixed. Completed work uses the OLD accepted coefficients and is never
repriced at the endpoint. This is first-order constitutive splitting; each
layer's explicit relative-viscosity-change bound requests time refinement.
Every source, coefficient, thermal and cancellation gate precedes the joint
publication, including when several normal/sliding islands share stores.
Existing thermal-network substeps remain the sole heat-transport owner.

Five focused groups in `tests/sliding_film_heat/feedback.rs` cover independent
two-mass Couette drag and layer work, distinct per-channel temperatures,
constant-source trajectory identity, mixed normal/sliding indices, source and
step-change refusals, complete heat-destination coverage and late cancellation.
The native coefficient-continuation regressions pass in the actual film-module
harness. Native execution of these coupled regressions remains unclaimed.

Four integration regressions cover independent two-body momentum and
discrete-gradient work, native reference equality and both Newton images,
the 0.8/0.2 Couette layer partition, partial selection and independent film
indices, original-network conduction, invalid route coverage/keys/duplicates,
and cancellation or late chart refusal with retry. Production and regression
sources pass Rustfmt parsing and bounded independent source review. Runtime
tests and a full crate build remain unexecuted in the full-storage environment.

