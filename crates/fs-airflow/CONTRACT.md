## Reynolds-dependent turbulent friction in native inertial graphs

The original `InertialGraph` accepts source-bearing `TurbulentPipeLoss`
branches alongside linear and quadratic losses. Each candidate midpoint uses
the original Colebrook pressure law's positive secant `Delta p(Q*) / Q*`.
The same port-Hamiltonian descriptor advances branch impulse and quadratic
kinetic energy; no startup friction factor, new momentum integrator or
constant-coefficient substitute is introduced.

The actual native midpoint is re-evaluated with the original nonlinear law
before the existing impulse and mechanical-energy gates can accept. The
reported average flow supplies transport, and the original pressure drop times
that flow supplies dissipated work. A material trial may refresh viscosity,
while retaining the original pipe geometry, roughness, law/source identity and
reference momentum density. Rechecking a thermal response preserves the
proposed native trajectory and tests its new physical residuals.

Initial, midpoint and final flows must stay in the authored Reynolds interval.
A sign reversal refuses because it crosses unsupported laminar/transition
flow, even if both endpoints are turbulent. Graph continuity, pressure ports,
resource limits, cancellation and immutable trial publication retain their
existing owners. A failed proposal or acceptance changes no accepted momentum.

Focused regressions cover signed mixed-network pressure/impulse/work against
an independent Colebrook solve, second-order unsteady decay and the steady
limit, viscosity refresh and stale-response rejection, and Reynolds/transition
refusal. This is developed, quasi-steady turbulent wall friction in an
incompressible fixed-reference-density graph. It does not supply turbulent
startup/transition, frequency-dependent wall-shear memory or compressible
acoustic storage.

# CONTRACT: fs-airflow

## Source-exact liquid channel energy

The material-bound liquid series, its solid-coupled solves and distributed
pressure-loss heating now use the original affine `Cp(T)` heat-exchange
integral. Their inlet conversion uses the original nonlinear `h(T)` inverse;
public enthalpy knots are not treated as an alternate piecewise-constant Cp
model. Compact source charts and optional finer boundary charts therefore
describe the same physical channel. Caloric traversal allocates no temporary
flow path and visits only original source intervals.

Hydraulic dissipation retains its declared constant-density, equal-area and
equal-elevation assumptions. The original `dp*Q` enters the fluid once, wall
heat remains distinct, and each Robin reference closes that same wall power.
Accepted enthalpy states still refuse increments lost at their absolute datum
or rounded across the driving wall temperature. Source support, source identity,
per-region balances and transactional refusal remain mandatory. Analytical
variable-Cp heating/cooling, non-knot inlets, datum shifts, region subdivision
and the original finite-solid power balance exercise the native consumer.

## Source-backed gas films for deposition

`graph::thermal::medium::mass_transfer::evaluate_laminar_plate_gas_film`
binds the original `BinaryDiffusionMaterial` to a thermal fluid owner. The
selected vapor must belong to its two distinct typed species; the other is
the declared carrier. Molecular diffusivity is resolved at explicit film
temperature and total pressure, preserving other source coordinates, selected
claims and receipts. Thermal and momentum views must agree at that temperature.
The generic fluid interface does not expose pressure or chemical identity, so
matching pressure/composition remains caller authority. Typed T/p aliases and
legacy untyped `T`/`p` coordinates refuse rather than retain stale source
applicability values. `at_film_state` refreshes the same selected diffusivity law;
source exhaustion cannot silently select another claim.

The source is Seader and Henley, *Separation Process Principles*, 2nd ed.,
Wiley, section 3.4, Eqs. (3-132), (3-137)--(3-138), available in the
[publisher's chapter](https://catalogimages.wiley.com/images/db/pdf/0471464805.ch3.pdf).
The model uses `Sc = nu/D_AB`, `Re = U L/nu`,
`Sh_L = 0.664 sqrt(Re) Sc^(1/3)` and `k_m = Sh_L D_AB/L`.
It admits positive speed, `Re <= 5e5`, and the conservative `Sc >= 0.6`
subset of the source's `Sc > 0.5` analogy. A flat, attached, thin laminar
boundary layer, zero streamwise pressure gradient, uniform interface
concentration and weak mass transfer are physical assumptions. Source support
does not establish those assumptions or propagate an uncertainty bound.

Surface strips retain their physical flow-facing leading edge. The stable
coefficient `k_[a,b] = k_[0,b]/(1 + sqrt(a/b))` avoids cumulative subtraction
and preserves integrated film transfer under subdivision at a common state
and concentration difference. Numerical cells do not restart the boundary
layer. Transport properties are frozen during each receiving phase step;
local variations require an independently admitted strip approximation.

`surface_kinetics` returns the existing thermal `SurfaceKinetics`, with
independently authored accommodation and this resolved film coefficient.
The reusable `mass_transfer::laminar_plate` owner also accepts explicitly
authored independent viscosity and diffusivity; the source adapter delegates
to that same checked geometry, correlation and kinetics handoff.
Original scalar coating, surface and chamber owners retain signed transfer,
finite inventory limits, latent/carried energy and recession. The caller binds
the same volatile species and compatible gas state to that owner. Its existing
gas-temperature pressure conversion and dilute linear film approximation are
unchanged. This bridge supplies no Lewis-number assumption, Stefan flow,
turbulent correlation, stagnant-gas diffusion geometry, nucleation or gas-flow
feedback. Tests cover actual source T/p, original claims and pair identity,
independent molecular transport, strip refinement, source/domain refusal,
cancellation and conserved deposition onto a bare substrate.

### Live original gas films in resistive and growing deposition histories

The mass_transfer::chamber adapter prepares a sparse, bounded roster of
SourceLaminarPlateChamberFilms. It pins each actual binary diffusivity selection
once, then re-evaluates that same D(T,p) law on every original full/half physical
trial and zero-time admission. Actual bulk total pressure is queried; volatile
partial pressure is not substituted. An explicit convention selects bulk
temperature or the arithmetic mean of the current deposit and bulk temperatures.
The current parcel's original model supplies viscosity/density at its evolving
composition. Bulk/equal-temperature evaluations reuse the parcel's existing
resolved state; a different film temperature evaluates that model at the bulk
pressure without heating or reconstructing a fictitious parcel.

The original material film-resolver seam drives fixed-conductance, arbitrary
enthalpy-chart conductivity and integrated conductivity-curve histories. It
retains the existing adaptive controller, molecular streams, separate condensed
heat stores, prescribed powers, gas reconstruction and acceptance accounts.
Film-bound growing layers must be the original planar shape: whole spherical
and cylindrical layers refuse. Strip positions remain measured from the same
authored physical leading edge; changing exposed area never invents a new
streamwise length.

A binary pair must contain the bound volatile and actual carrier; the carrier
must remain positive and foreign active gas species refuse. A carrier-only
initial state is accepted without adding a vapor seed. Only canonical T/p
coordinates are allowed in the binary source context, so a composition or other
state coordinate cannot silently remain frozen. Each binding declares a maximum
volatile mole fraction in (0,0.05]. Every source evaluation checks the actual
bulk mole fraction, saturation pressure divided by total pressure, and the
original two-temperature kinetic-equilibrium pressure divided by total pressure.
This five-percent cap is a conservative implementation admission policy, not a
published universal accuracy bound or proof of weak blowing/Stefan flow.

The prepared source owner is immutable during trial evolution. Coefficients
reuse the material owner's scratch map, and each bound key is overwritten on
every invocation. Source/domain/correlation/dilute-limit failures preserve their
typed error and abort without a physical candidate. Cancellation is polled per
bounded row and between original constitutive queries. The existing continuum,
attached laminar, dilute weak-transfer and frozen-strip approximations remain;
no turbulence, Stefan/Soret correction, nucleation, moving velocity field or
additional gas-side heat law is inferred.

### Live original gas films in liquid evaporation

`mass_transfer::evaporation::SourceLaminarPlateEvaporationFilms` binds those
same molecular-diffusivity, mixture-transport and plate/strip owners to the
original simultaneous liquid/thermal/finite-gas endpoint and adaptive event
timeline. Each selected native cell retains its liquid law, volatile identity,
caloric binding, admitted enthalpy interval and sensible heat coefficient.
The prepared binary source is immutable; every actual nonlinear and adaptive
trial refreshes its selected diffusivity at current total pressure and the
explicit bulk or mean liquid/gas film temperature. The actual current mixture
supplies viscosity and density. An original source-domain exit refuses even
when another source claim could cover the trial.

Only mass-transfer coefficient and representative film temperature are
replaced. The coupled owner retains its exact total and volatile partial
pressures and validates the complete active surface roster after each callback.
Its Stefan/interface law, mass balance, latent energy, upstream vapor enthalpy,
reciprocal sensible heat and original gas caloric reconstruction remain the
only physical owners. The source adapter adds no heat or vapor seed. Source
failures preserve their typed cause and publish no thermal or gas endpoint.
The adaptive timeline re-evaluates the same sources in coarse/fine/rejected
trials, accepts only fine physical transfers and clips at native events.

Bindings may be sparse; unbound interfaces retain their authored films. Bound
cells must be active at preparation and remain present. Unbound births retain
the original timeline semantics. The source context admits only canonical T/p
axes; foreign active gas species, zero carrier, nonplanar declarations, strong
vapor loading, stagnant flow and correlation-domain exits are unsupported.
Planarity, attached thin flow, continuum validity and weak normal transfer
remain caller declarations. Both bulk vapor fraction and saturation pressure
fraction must satisfy the existing authored ceiling in (0,0.05]; this policy
does not validate a high-blowing boundary layer. There is no new droplet,
turbulent, plume, gas-velocity or resolved liquid-interface-temperature model.

Focused G0/G1 tests exercise live original pressure/composition/temperature
queries, retained selected claims, both signed phase directions, a direct
Stefan-rate reconstruction, independent gas/liquid energy and species balances,
wind-speed response from a carrier-only atmosphere, adaptive heater events,
source exhaustion without claim switching, cancellation and atomic refusal.

### Original molecular sources around spherical liquid droplets

`mass_transfer::evaporation::droplet::SourceSphericalDropletFilm` pins one
explicit binary molecular-diffusivity claim and volatile species. Every phase
or heat/phase interval samples that same claim at the actual initial bulk gas
temperature and total pressure. The original spherical geometry supplies D/r;
no viscosity, thermal-diffusivity analogy, forced-convection correlation or
flat-plate dilute cutoff is introduced. Current active gas constituents must
belong to the authored pair, a represented positive carrier must remain, and
each call's liquid binding must retain the pinned volatile identity. An absent
initial vapor needs no artificial seed. Molecular source contexts admit only
canonical temperature and pressure coordinates.

The adapter delegates original droplet phase and reciprocal sensible-heat
owners, retaining IF97 or other original liquid calorics/density, NASA gas
calorics, selected interface kinetics, Stefan flow, finite inventory and all
mass/energy accounts. The original first-order interval freezes initial area
and D/r, including the heat-then-phase split. Actual final gas and size are
qualified against the same selected claim before publication; a source exit
refuses even if another claim covers the endpoint. Initial/final source-film
samples expose those actual coordinates and receipts. The refreshed final D/r
is a next-interval diagnostic, not a second transfer or an implicit correction.
Geometry remains allocation-free; material-source resolution retains its
existing receipt allocation behavior, without an additional film map or roster.
Cancellation and source/phase refusal publish no partial physical candidate.

Added tests cover actual T/p and size response, original source identity,
carrier-only and nondilute binary states, volatile/context/source refusals,
original phase/heat candidate equivalence, next-interval refresh and late
same-claim exhaustion. The original IF97/NASA integration uses unchanged
checked-in NASA coefficients and explicitly manufactured diffusivity data;
it establishes conservative interoperability, not experimental drying accuracy.

## Original coolant inlet, wall, flow and conductance linearization

`CaloricAirPath::linearize` prepares the complete inlet/wall temperature
response during one native source-Cp march. It retains each original
segment's enthalpy gain, stable gain complement, wall enthalpy response and
conductance-weighted mean response. All caloric derivatives use the original
variable Cp at the actual fluid states; wall reservoirs need no extrapolated
material law. The prepared operator borrows the immutable path and wall
operating point, and refuses actual fluid endpoints on the source boundary
where a two-sided derivative is not admitted.

`CaloricPathLinearization::apply` propagates one inlet/wall direction to every
outlet temperature, heat rate and Robin reference, plus the final outlet and
total heat. `apply_transpose` returns inlet and ordered wall gradients for
arbitrary seeds on those same outputs. The two recurrences take O(segments)
work, use caller-owned output buffers and allocate nothing on success. They
perform no material re-query, caloric inversion, finite difference, dense
Jacobian assembly or additional nonlinear solve. Preparation owns O(segments)
coefficients and the returned native march.

`apply_with_flow` additionally accepts a prescribed mass-flow direction in
kg/s and one h*A direction in W/K per original segment.
`apply_transpose_with_flow` returns their objective gradients along with the
original inlet/wall gradients. Both differentiate the same source-integrated
enthalpy equation: the flow response includes the direct change of transported
mass as well as residence-time feedback. The stable original moments retain
`T_wall-T_out` and `T_out-T_mean` before absolute-temperature rounding; weak
flow sensitivity and saturated conductance sensitivity are not recovered by
subtracting nearly equal heat transfers. These actions retain linear work,
caller-owned scratch outputs and allocation-free repeated application.

Weak heat derivatives use the retained gain complement instead of subtracting
nearly equal endpoint enthalpy tangents. Total heat tangents use a compensated
sum of the same local watt changes. The original positive mean-response
moments retain small wall sensitivity and the exact equilibrium limit.
Cancellation and invalid input return no physical state change; action output
is explicitly scratch storage and may be partly written on error.

This is the derivative of the original single-phase channel at fixed pressure,
composition and source law. Flow and h*A are independent prescribed inputs;
hydraulic, geometric or correlation owners must supply any design-to-input
derivatives. Solid response, return-loop controls and coupled design adjoints
retain their own owners; they are not inferred by this channel operator.
Regressions cover independent variable-Cp channel integrals, inlet and wall
actions, transpose gradients, weak transfer, equilibrium, source admission
and cancellation. Flow-control regressions additionally compare independent
variable-Cp channel integrals, every control gradient, the Robin heat tangent,
simultaneous flow/conductance scaling and weak/equilibrium/saturated limits.
The original numerical kernel and all five action regressions pass in a focused
harness with fixture path/material containers and execution checkpoints;
whole-crate integration and the native cancellation runtime remain unverified.

## Source-caloric thermal networks

`graph::thermal::transport::caloric::CaloricTransportNetwork` binds one
immutable `IsobaricHeatCapacityLaw` to the original solved hydraulic DAG or
closed, thermally anchored hydraulic circulation.
All supplies, branch exchanges and junctions share the same fixed-composition
isobar. Junctions mix mass-weighted relative enthalpy and recover temperature
with the original source inverse. Branches reuse `CaloricAirPath`, including
its retained weak-exchange heat and directly integrated Robin mean.

The original graph walk, actual flow orientation, region ownership, cancellation
and whole-network watt ledger remain shared with constant-Cp transport. The
ledger advects retained enthalpy, including hydraulic mass residual energy;
it does not reconstruct branch heat from rounded outlet temperatures. The
explicit physical reference temperature selects the Cp accompanying the
fixed density snapshot independently of the source's arbitrary energy datum.

`solve`, `solve_from`, `solve_iqn` and `solve_iqn_from` use the existing
shared-solid driver and source-aware branch balance probes. Fixed/Aitken/IQN
acceleration, temperature convergence and per-region watt gates retain their
native owner. Source/domain or balance refusal returns no partial field.

`with_recirculation` admits the native bounded, prescribed return-link topology
with receiving-supply mass fractions and strictly positive fresh makeup.
It solves the nonlinear mixer equations in relative specific enthalpy, using
the original channel inlet response and native mass-weighted DAG elimination.
The reduced Newton solve has explicit iteration and damped-trial limits;
every trial remains inside the original source support. The original shared
solid driver consumes the resulting references without a second coupling loop.

Supply enthalpies are primary unknowns. The shared graph walk retains inlet
and node enthalpy coordinates, so returning exhaust is never reconstructed from
rounded temperatures. Mixer watt residuals and the source-minimum-Cp bound on
temperature residual must both pass. A separate fresh/undrawn-exhaust ledger
then closes the outer heat balance. Reports retain actual mass rates and mark
the physical reference temperature of diagnostic Cp capacity rates; those
rates do not enter the nonlinear heat balance.

`new_circulating` admits closed native hydraulic loops and simultaneous
mass-weighted node enthalpy balances. Pressure boundaries are hydraulic datums,
not external fluid supplies. Original direction, continuity, reference-state,
64-node and resolved thermal-anchor checks remain with native admission.
Every nonlinear trial recomputes source channel heat and analytic inlet gains.
A bounded damped Newton solve requires both mixer watt closure and a small
correction measured through the minimum source Cp; a small heat residual alone
does not establish convergence under weak wall exchange. This is a local
convergence check, not a certified nonlinear or physical error bound.

Closed-loop residuals retain compensated relative enthalpy differences and
original branch heat. The final native march independently repeats mixer and
whole-network watt admission. The supported wall-temperature midpoint is an
initial guess; an unsupported trial may fall back once to the admitted physical
reference. Neither guess is a reservoir. Bath temperatures may exceed source
support only while every fluid state and trajectory remains supported. The
existing shared-solid drivers consume the same accepted loop response.
Imposed return links cannot be added to a closed circulation.

Hydraulic flow, density and film coefficients remain prescribed. This path is
steady and single-phase. Composition mixing, return-path fan/pressure work,
density/flow feedback and their adjoints are not supplied by this specialization.

Focused regressions exercise real hydraulic split/merge and reverse-flow
solutions, nonlinear mixing, source/datum invariance, weak heat, shared-solid
power balance, hydraulic defects, cancellation and source refusals.
Return regressions include an independent affine-Cp analytical balance,
heating/cooling, 99.99% recirculation, multi-supply mixing, source-only Cp
changes, datum/order invariance, bounded refusal, a coupled solid and retained
weak heat at a 1e-14 W balance tolerance. These are numerical model checks,
not experimental validation or a global nonlinear-convergence guarantee.

## Source-specific heat in steady conjugate channels

`conjugate::caloric::CaloricAirPath` binds the original
`IsobaricHeatCapacityLaw` to native ordered `AirSegment` geometry and a
prescribed mass flow. It integrates `Cp(T)/(T_wall-T)` over the traversed
source intervals to obtain each outlet, then reports heat from the original
relative enthalpy integral. Source knots, pressure/composition context,
identity, receipts and uncertainty remain with the immutable material law.
Cloning shares that law and its prepared caloric index.

The relative temperature increment is retained through the original heat
integral, and the conductance-weighted mean temperature is integrated directly.
This keeps weak exchange accurate without amplifying rounded absolute-outlet
temperatures through division by a tiny conductance. Nonzero exchange below
the representable outlet-temperature resolution still refuses.

The Robin reference uses the same actual watts as the stream. `solve`,
`solve_from` and `solve_iqn` reuse the existing conjugate driver's response
admission, region ownership, fixed/Aitken/IQN acceleration, cancellation,
temperature tolerance and independent per-region watt gates. Original
constant-Cp goal and refinement formulas cannot silently consume this distinct
path type. Per-segment NTU is the inlet-Cp diagnostic; effectiveness comes from
the integrated temperature change, not an assumed constant-Cp exponential.

`march_with_inlet_response` returns the analytic fixed-wall inlet response of
that same march. The original material exchange retains `dh_out/dh_in` before
rounding absolute temperatures, including the exponential equilibrium limit
and saturated outlets. Channel composition multiplies these passive enthalpy
gains; temperature response additionally uses the original inlet/outlet Cp
ratio. This is a derivative at fixed mass flow, conductance and wall states,
not a coupled design adjoint. Analytical constant/affine-Cp tests, central
differences, heating/cooling, source knots and datum/refinement checks cover
the response.

The model remains steady, single-phase, one-dimensional and isobaric, at fixed
composition, mass flow and heat-transfer coefficient. It adds no hydraulic
feedback, pressure work, temperature-dependent h, mixing junctions or a coupled
solid/hydraulic adjoint.
Walls may be outside the fluid's caloric support; the actual traversed fluid
temperatures may not. A source refusal returns no partial march or solution.
Tests compare constant/affine-Cp analytic channels, weak exchange, heating and
cooling across source knots, segment refinement, datum independence, support refusal and a
power-driven coupled solid/stream balance.

## Transient original Carreau hydraulic networks and storage

`graph::rheology::transient::TransientCarreauGraph` connects original
`RadialCarreauPipe` histories through ideal pressure junctions and explicitly
authored linear hydraulic storage. Every branch retains its complete annular
velocity field, original thermal Carreau law and rate ceiling. Density is fixed
and each branch's prescribed temperature is held during each interval. The
developed section remains isothermal; different branches may have different
temperatures. This consumer adds coupled dynamics
to the existing Carreau topology; it does not fit an effective resistance or
replace the radial field with a mean-flow inertance.

All branches advance from the same accepted initial histories and duration.
Unknown midpoint pressures satisfy `Q_out - injection + C*(p_next-p_old)/dt = 0`,
with `p_next = 2*p_mid-p_old` at storage nodes. The positive branch derivative
comes from the accepted original radial tangent. The existing checked-CG
conductance Laplacian gains only the positive `2*C/dt` storage diagonal, with
the same Jacobi preconditioner and checked residual. Bounded Newton and line
search retain every original source bound. A separate branch-solve budget
counts all attempted original solves, including rejected pressure trials.

An absent storage law denotes a zero-volume ideal junction. Its initial
radial flows must satisfy the held injections, and accepted midpoint and final
flows must both meet continuity tolerances. A discontinuous prescribed-flow
change at such a node refuses rather than inventing an impulsive pressure.
Each connected component retains the existing explicit pressure-reservoir
anchor requirement. A storage node cannot also be an ideal pressure reservoir.
Each positive constant compliance has a declared pressure interval; both
physical storage endpoints must remain inside it. Storage pressure uses the
same reference as all pressure ports, with incremental energy `C*p*p/2`.

Reports retain original per-branch kinetic-energy change and nonnegative
viscous heat, storage volume and energy change, signed reservoir work and
prescribed-flow work. Free-node residual work, original radial work defects
and midpoint storage roundoff remain separate from heat. An arithmetic work
identity is checked before publication; its tolerance is not a source,
discretization or physical-error certificate. Trial radial histories and
evaluation buffers are reusable. Existing CG corrections and owned returned
reports still allocate graph-sized storage. Only a complete globally admitted
step swaps all branch histories, storage pressures and common time. Source
refusal, work exhaustion, invalid storage and cancellation preserve accepted
history and permit an identical retry.

`TransientCarreauGraph::at_temperature` continues this same network at one
prescribed uniform temperature between intervals; `at_temperatures` accepts
one local temperature per original branch in branch order. Every original pipe must
admit its unchanged physical rate ceiling against the shifted reference support;
a refusal in any branch or cancellation returns no partially rebound graph.
The original topology, density, annular masses, compliance laws and compatibility
identity remain unchanged, so accepted radial histories, storage pressures, time
and existing workspaces carry forward. Preparation copies graph metadata but
neither radial history nor solve buffers. The actual new nonlinear viscosity
enters each subsequent original pressure trial and physical heat calculation.
This requires the authored fixed-density and fixed-compliance approximations
to remain valid; temperature evolution, thermal energy and expansion need their
own owner. No thermal mixing or axial temperature profile is inferred.
G1/G3/G4 tests exercise retained storage/radial history, both old and new
workspace compatibility, original uniform and regional source trajectories,
an independent dense radial/compliance matrix with local Arrhenius viscosities,
cancellation and atomic refusal of a later branch's reduced-rate ceiling.

Each transient connected component may instead be sealed and anchored by its
own positive physical storage. Stored pressures remain unknowns in the original
Schur solve; no storage node becomes a pressure boundary and no artificial
reservoir, fixed gauge or volume sink is introduced. Net signed injection changes
the actual stored volume, while zero-storage junctions retain their original
continuity constraint. Disconnected sealed and reservoir-driven components may
coexist. A component with neither storage nor prescribed pressure still refuses,
and steady Carreau solves still require a reservoir in every component.

`TransientCarreauGraph::step_with_heads` adds one explicit held pressure rise
per original directed branch. The same resolved radial law receives
`p_from-p_to+head`, so a source can drive circulation in a sealed compliant loop.
Its signed work is `head*Q_mid*duration` on the accepted trajectory, reported
separately from reservoir/injection work, viscous heat and kinetic/compliance
storage. Opposing heads can recover energy from retained fluid inertia. All
loads participate in the existing coupled pressure solve and original energy
accounting; this supplies no inferred pump curve, shaft efficiency or cavitation
law. The ordinary `step` keeps its original zero-head arithmetic path.

The model is a radially resolved, axially developed, laminar, fixed-density
hydraulic network with lumped incremental compliance. It does not infer a bulk
modulus or nominal fluid volume, resolve axial acoustic waves/water hammer,
model junction/entrance losses, advance thermal state or predict cavitation.
Implicit midpoint still requires temporal and radial resolution; it can
oscillate when stiff modes are under-resolved. G1/G3 tests compare a coupled
Newtonian velocity/storage step against an independent simultaneous dense
system, preserve the original nonlinear series trajectory, resolve loop
continuity and reversal, differentiate the nonlinear Schur tangent, retain
different radial histories at equal mean flow, and check separate storage,
injection and reservoir work. Source, partial-work and cancellation refusals
leave the complete state unchanged. Existing steady Carreau regressions apply
to the shared pressure operator.

## Self-heating of retained Carreau momentum

`graph::rheology::thermal::transient::TransientThermalCarreauGraph`
couples the original radial Carreau histories and hydraulic compliance to one
explicitly shared, well-mixed `TabulatedHeatStorage` inventory. Its capacity
includes all pipes, compliant volumes and returning reservoir fluid. Every
thermal trial starts from the same accepted velocities, stored pressures and
clock, prepares original viscosity at the trial midpoint temperature, and
performs the original coupled transient pressure/radial solve.

Only the native trajectory's integrated nonnegative viscous heat enters
sensible storage. Pressure/reservoir/injection/head work, radial kinetic energy
and hydraulic compliance energy retain their separate original owners. Thus
startup can store pump work as momentum, coastdown can heat the fluid without
pressure input, and an opposing pressure head can recover mechanical energy
while the fluid still dissipates heat. The existing bounded scalar caloric
midpoint bracket/refinement is shared with the quasisteady thermal owner;
original source support and physical rate ceilings remain unchanged.

`advance` and `advance_with_heads` borrow the accepted history and return the
complete next temperature owner, radial/storage state and work report together.
The same original workspace remains compatible across accepted temperatures.
Failed trials, endpoint source admission, energy refusal, exhausted shared
hydraulic/branch work and cancellation publish no partial physical state.
Both sensible-energy and combined mechanical-plus-thermal residuals must meet
declared tolerances; those numerical gates do not certify discretization error.

Temperature is uniform within this entire declared inventory; there is no
spatial enthalpy advection, finite mixing time, open-stream donor temperature,
density change, thermal expansion, or inferred fluid mass/heat capacity.
The scalar continuation selects a bounded local heating/cooling branch and
does not certify uniqueness or unresolved equilibria. Radial and temporal
refinement remain necessary. Independent two-annulus momentum/caloric equations
cover startup, coastdown and recovering heads; a simultaneous pressure/storage
oracle separates compliance energy from heat. Zero-flow cooling, the balanced
steady limit, original-workspace continuation and atomic retry are also tested.

## Regional enthalpy transport through original Carreau networks

With `transient-transport`,
`graph::thermal::transport::transient::carreau::RegionalCarreauTransport`
couples native regional thermal storage and conservative enthalpy transport
to the original quasi-steady Carreau pressure network. Each real pipe section
has one distinct native fluid cell; each hydraulic junction has additional,
explicit finite mixing inventory. Pipe length and volume must match the original
centreline length and `pi R^2 L`, and incident pipe endpoints agree at junctions.
All fluid masses equal the declared fixed reference density times their volumes.
The original registered enthalpy material, declared caloric identity, complete
geometry and masses remain pinned. Native wall/solid cells may coexist outside
the fluid map; no heat capacity, material source or junction volume is inferred.

Every native backward-Euler trial resolves each pipe's OWN current cell
temperature through the original caloric law, prepares its unchanged original
Carreau-Yasuda source and physical rate ceiling, and solves the original network.
Its actual signed volume rates drive native junction-to-pipe-to-junction mass
transport. Native donor selection reverses with flow and retains heated fluid
already stored in each pipe. External streams attach explicitly to junctions;
pressure datums supply no implicit mass or donor enthalpy. The native roundoff
mass-closure gate remains stricter than an arbitrarily loose hydraulic tolerance.

Each original accepted steady profile's nonnegative viscous `dp*Q` is deposited
once in its own pipe cell. Reservoir/injection pressure work and constitutive
or continuity residual work remain in the original hydraulic report and are
never added again as heat. Internal enthalpy transfer stays internal; native
net external stream enthalpy is reported separately. Caller-owned wall links,
ambient coefficients and additional sources enter that same original residual.
The existing nonlinear implicit driver is shared with prescribed-flow transport;
it owns bounded endpoint evaluations, native caloric inversion, cancellation,
independent final evaluation and local/global energy acceptance. Refused source,
inventory, stream, work or cancellation proposals publish no changed state.

This is transient regional heat with quasi-steady developed laminar momentum.
Each authored section is thermally well mixed, and its endpoint flow represents
the full implicit interval. Time and spatial refinement remain required;
momentum history, radial thermal gradients, entrance effects, buoyancy, expansion
and the instant of flow reversal are outside this owner. The immutable graph
accessor is the original template; accepted profiles and temperatures are in
the returned endpoint report. The retained-radial-history thermal owner remains
the separate explicit one-inventory model described above.

Six focused regressions cover independent serial Poiseuille/Arrhenius and
finite-inventory energy roots, reverse flushing with nonlinear heat capacity
and a large enthalpy datum, original nonlinear Carreau feedback, the zero-flow
heat limit and parity with prescribed-flow stepping, strict admission/work
refusals and cancellation after an actual hydraulic solve. All 54 tests in the
source-linked Carreau/thermal battery passed with original material, radial,
advection, backward-Euler, checked-CG and Newton/Krylov numerical code. The
execution context and unrelated module namespaces were isolated to avoid the
unavailable full dependency build; this is not a full-crate/runtime test pass.

## Transient radial Carreau momentum with regional enthalpy transport

With `transient-transport`,
`graph::thermal::transport::transient::carreau::transient::TransientRegionalCarreauTransport`
composes the original retained radial `TransientCarreauGraph` with the regional
finite-inventory binding above. The constructor requires the same original
pipe topology, source laws, rate ceilings, geometry and reference density.
Every hydraulic compliance is explicitly refused: fixed thermal cell masses
cannot absorb pressure-dependent storage volume. The separate one-inventory
thermal owner remains the compliant-network path.

Every native backward-Euler thermal trial starts from the SAME accepted radial
histories. Each original pipe is prepared at its own candidate ENDPOINT
temperature, then its original midpoint momentum step resolves pressures,
annular velocities and integrated irreversible heat. This is endpoint-held
rheology with midpoint momentum, first order in coupled time; no midpoint
temperature is substituted or claimed. Native midpoint quadrature of mean
volume flow drives conservative signed junction/pipe enthalpy transport.
Only original integrated viscous heat divided by the interval duration enters
each pipe cell's heat source. Steady `dp*Q`, pressure work, kinetic-energy
change and hydraulic numerical residuals are never substituted for that heat.

Explicit external streams use the actual native mean rates. Pure callbacks
receive original old, mean and final flow diagnostics to admit source-dependent
exchange. The original native external-power routine evaluates other authored,
ambient and reservoir heat before advection or viscous heat is added.
Independent thermal and mechanical gates precede the combined
`delta(thermal + kinetic) - reservoir - injection - head work - stream enthalpy - other heat`
check. Reports keep all those physical terms, native numerical residuals and
aggregate branch-solve counts separate. Initial physical clocks must agree;
the next prepared owner, complete radial state and native thermal timeline
are returned together. Source, inventory, stream, clock, work and cancellation
refusals preserve both accepted physical inputs; reusable workspaces remain
scratch and can be retried.

Each developed section remains thermally mixed and fixed in geometry/density.
The native mean is a midpoint quadrature, not an exact continuous trajectory
integral. Within-interval reversal times and opposing traversals are unresolved;
time refinement is required for their residence and heat transport. The original
developed radial model still omits axial momentum/kinetic-energy convection,
junction kinetic losses, thermal expansion and radial thermal gradients.
No compressible or changing-mass interpretation is implied.

Six new regression sources cover an independent two-annulus momentum and
regional caloric root during startup/coastdown/recovering head work, distinct
nonlinear regional source temperatures, zero/steady limits, equal initial mean
flow with different retained radial heating, refusal/retry and cancellation.
These six regression cases now pass in the original-source slice described
below; three fixture calls were corrected to pass the native owned velocity
roster. The earlier 54-test result applies to the preceding quasi-steady delivery.

### Resolved regional transport through reversal

`carreau::transient::adaptive::advance_resolved` resolves an event-free window
by comparing one original coupled interval with two native half intervals.
Only independently admitted fine states and reports are retained. Every native
thermal cell and every radial annulus participates in the local comparison;
mean flow alone cannot hide a different shear history. Relative thermal
tolerances scale interval energy changes rather than the formation-energy datum.

Separate limits bound the duration of sign-changing native intervals and
`min(V+, V-)` of each branch's affine numerical flow trajectory. A near-zero
mean therefore cannot bypass requested counterflow resolution even at uniform
temperature. Chronological original substeps retain both directional transport,
local viscous heat, kinetic change, external streams and signed pressure/head
work; discarded comparisons contribute only to the aggregate work count.
No new caloric integrator, effective viscosity or extrapolated physical state
is introduced.

Native attempts, retained steps and total branch solves have explicit limits.
Original per-step thermal/nonlinear/source gates remain active; each original
branch-work cap is clipped to the remaining window allowance. Source, mass,
numerical and conservation refusals are terminal, including discarded trials.
The complete window publishes synchronized owners together or preserves both
input histories, even after a locally accepted prefix; scratch remains reusable.

This is first-order local error estimation and resolution of sampled affine
reversals, not a rigorous global bound or exact event locator. Held mechanical
ports, original fixed mass/density and spatial applicability remain unchanged.
Split forcing discontinuities and native thermal events explicitly, and refine
maximum interval and spatial inventory for the intended accuracy.

The focused original-source slice passed 53 tests, including the six regional
transient cases and five new adaptive cases: independent coupled annular/caloric
ODE refinement, uniform-temperature counterflow, distinct equal-mean radial
histories, complete-window work/source refusal with retry, and cancellation
after a local accepted prefix. It executes original material, radial, thermal,
advection and nonlinear solver code with context-only cancellation glue; this
is not a full Cargo or native-runtime verification. Pinned rustfmt and focused
diff checks passed.

### Authored pressure, head and injection histories

`carreau::transient::adaptive::advance_resolved_with_drive` uses the same
native adaptive controller as the held-drive entry. A pure callback supplies
original pressure boundaries, volume injections and optional branch heads at
each proposed interval's strict interior midpoint. Each coarse and half-step
proposal gets its own single query; every thermal residual and derivative probe
within that native step reuses the same held forcing. No waveform, interpolation
or substitute pump model is inferred.

Finite strictly increasing authored discontinuity times clip comparison windows
exactly. Neither coarse nor fine proposals cross an event, and no callback is
queried at the discontinuity itself. The next forcing segment starts there with
the actual retained radial and caloric histories and a fresh comparison-width
proposal. Nonbinary event periods do not carry a rounded previous period into
a final unrepresentable time sliver. Native thermal inventory events retain
their separate admission.

Accepted `drive_samples` correspond one-for-one with chronological native
reports and retain their start, end, midpoint and exact physical port values.
The exchange callback receives that same sample to author actual incoming and
outgoing streams, including separate injection/reservoir streams at one node.
Rejected comparison samples carry no physical work or heat credit. Drive calls
are counted and bounded by native attempt limits; source refusal, cancellation
and exhausted work preserve the complete input window even after a locally
accepted prefix.

Pressure and head waveforms, and changing injections at pressure-reservoir
nodes, retain original midpoint work accounting. A free rigid junction still
requires the native held injection to match initial, mean and final flow.
Incompatible changes refuse even at an event; no impulsive pressure, velocity
reset or changing-mass storage is fabricated. The coupled thermal method remains
first order, with explicit maximum-interval and source-support limits.

Five new regressions check exact held-drive parity, nonbinary event splitting
against separately executed native segments, smooth pressure/head/injection
forcing against independent radial/caloric dynamics and mechanical work,
invalid/late drive and rigid-junction refusals, retry, cancellation and bounded
queries. The same original-source slice now passes 58 tests; full Cargo and
native-runtime validation are not claimed. Pinned rustfmt and focused diff
checks passed.

## Original-law Carreau network design derivatives

`CarreauPipeGraph::solve_linearized` and
`solve_with_injections_linearized` retain a successful original radial solve
inside an immutable `CarreauGraphLinearization`. No external or mutable
`CarreauGraphSolution` can admit cached conductances. Forward directions and
scalar-objective adjoints cover each pipe's log radius, log length and prescribed
temperature, plus original-order pressure boundaries and signed flow injections.

With accepted signed flow Q, constitutive pressure drop p and original positive
conductance G=dQ/dp, radial momentum gives dQ/dln(R)=3Q+pG and dQ/dln(L)=-pG.
For the source's original logarithmic shift slope s=dln(a_T)/dT, the Vertical
calibration gives dQ/dT=-s pG; VerticalAndRate gives dQ/dT=-s Q. These distinct
responses retain the authored temperature convention. No fitted apparent
viscosity, radial finite difference or derivative of Newton iteration history
is substituted. Source-law and physical rate ceilings must contain an open
neighborhood, including the complete prepared rate envelope under a horizontal
temperature shift; a boundary state refuses a two-sided derivative.

The same grounded conductance Laplacian and checked CG used by the primal solve
each direction or adjoint. Original profiles, parameter responses, Jacobi
diagonals and component references are retained and reused. The coefficient
storage and each operator application are O(nodes+pipes); one adjoint handles
all exposed parameter gradients. Rest uses the true finite conductance.
Pressure directions are centered on one retained boundary per component before
the linear solve, keeping large uniform datum changes out of the shear-flow
accuracy scale. Pressures are then restored and flows retain the centered
calculation.

Differentiated continuity is actual pipe outflow minus authored injection.
Repeated inputs share their node derivative while retaining input order;
injection at prescribed pressure changes net reservoir exchange without changing
the hydraulic state. Each tangent reports the physical continuity residual and
the checked linear residual. Adjoint weights act on node pressures and branch
flows; direct objective dependence on geometry, temperature or ports remains
the caller's chain rule.

Focused regression sources cover independently differentiated
Arrhenius/Poiseuille flow, original solved-network differences for radius,
length and temperature under both shift modes, forward/adjoint dot products,
repeated and fixed-node injections, stagnation, component gauges, material
support boundaries, exhausted work and cancellation/retry. The changed Rust
files were parsed/formatted; native Cargo compilation and tests were not run
in the unavailable constellation workspace. Independent numerical algebra
checks are not a native-runtime or performance campaign.

These are local nominal-model derivatives. The primal constitutive/quadrature
error and source uncertainty remain, and a checked linear residual is not a
gradient-error certificate. No finite perturbation, changed laminar regime,
transient momentum, self-heating evolution or coupled thermal adjoint is
admitted by this derivative owner.

## Steady shear-thinning hydraulic networks

`graph::rheology::CarreauPipeGraph` connects original `LaminarCarreauPipe`
owners through arbitrary ideal pressure junctions, including parallel paths,
loops, flow reversal and disconnected components with their own reservoirs.
Every branch retains the same original thermal Carreau-Yasuda law; geometry,
prescribed local temperature and explicit laminar wall-rate ceilings may differ.
Each developed pipe remains isothermal within its own section. Volume continuity
assumes one incompressible fluid at fixed density, with no variable-density
mass-flow, junction heat mixing or thermal evolution inferred.

`CarreauPipeGraph::at_temperatures` prepares all original branch laws at the
authored local temperatures without replacing a source law, geometry or physical
rate ceiling. Malformed rosters, source refusal or cancellation publish no
partial owner. The resolved local viscosities and differential conductances
change the actual pressure/flow split and viscous heat. An independent
Arrhenius/Poiseuille loop fixture checks pressures, branch flows, tangents,
source retention and power. The single-inventory `ThermalCarreauGraph` explicitly
requires equal branch temperatures; it refuses a nonuniform hydraulic graph
rather than discarding regional thermal state.

The regional-temperature and transient caloric changes passed the pinned
formatter, whitespace checks and a 48-test source-linked Carreau battery.
That check compiled the original material, radial, graph, checked-CG and caloric
algorithms and their real regression sources. Execution-context glue was isolated
because the native crate dependency build could not complete; no full-crate or
asynchronous-runtime integration pass is claimed.

`solve_with_injections` also admits signed prescribed nodal volume ports:
positive inputs supply fluid and negative inputs withdraw it. The required
free-node pressures follow from actual pipe outflow minus authored injection,
using the same radial material response. Repeated inputs at a node retain
their original list and are combined by deterministic compensated summation.
A pressure reservoir may coexist with an injection; its net exchange is the
actual pipe outflow minus that injection. Such a port does not change the
pressure-controlled pipe solution or loosen its convergence tolerance.

Each pressure trial resolves the original radial viscosity profile. Junction
Newton corrections use that accepted profile's positive differential
`dQ/d(delta_p)`, including its finite true zero-flow limit, in a grounded
matrix-free graph Laplacian. The existing checked CG and Jacobi preconditioner
solve each correction; no apparent-viscosity replacement, secant Jacobian,
epsilon resistance or independent Krylov implementation is introduced.
Storage is O(nodes + pipes). Accepted profiles and conductances are reused
through the linear solve and retained unchanged in the final result.

Bounded pressure/injection continuation and damped corrections can recover when an
intermediate guess exceeds a pipe's source support even though the final
network is admissible. Each anchored component starts at its own uniform
pressure, so another component's pressure gauge cannot consume its rate
support or continuation work. Every evaluated pipe stays under its original source
and numerical controls. Demand magnitudes enter free-node continuity scales.
Only convergence at every original target pressure and injection may publish;
unresolved source states, exhausted work, unreferenced
components and cancellation return no partial network. Checkpoints surround
each finite pipe evaluation and bounded graph/linear traversals.

The returned signed net reservoir work, separate prescribed-flow pressure work,
nonnegative accepted viscous dissipation,
constitutive pressure-residual work and remaining free-node residual work
satisfy an explicitly reported floating-point power balance. Numerical defects
are never relabeled as heat or silently omitted. Heat is available to a thermal
owner and is not deposited by this steady solver. Continuity convergence does
not certify radial quadrature accuracy, source uncertainty or a laminar regime.
The original no-slip, fully developed pipe assumptions remain; there is no
new gravity, transient inertia, entrance, turbulence or compressibility model.
Balanced injection/withdrawal pairs represent ideal flow-controlled transfers;
their pressure work can be supplied or absorbed. No pump head curve, efficiency,
shaft state or transient source dynamics is inferred.

## Unsteady laminar pipe networks with internal continuity

`graph::inertial::laminar::network::LaminarNetwork` connects original
`LaminarPipe` components through arbitrary ideal internal junctions,
branches and loops. A boundary-rooted spanning forest supplies integer
through-flow/circulation columns. Orthonormalizing their original inertia
metric eliminates only redundant mean-flow coordinates; every radial
memory arm remains. Internal volume continuity belongs to the continuous
state space, not an endpoint or interval-average repair. Unresolved
numerical rank refuses instead of discarding a physical circulation mode.

The projected system retains original positive viscous source squares.
`PreparedExactZoh` advances one coupled held-input trajectory, and
`PreparedSourceLoss` measures each branch's actual heat on that trajectory.
Reservoir ports and explicit pressure-head ports remain power conjugate to
their original flows; the latter also provide interval-average branch flow
without duplicate observer exponentials. Reservoir and head work are
reported separately. Gauges are normalized per physical component before
matrix products. Work and forcing roundoff bounds use uncancelled port
terms so recovering and driving sources may coexist.

Original fixed-density sources, branch Reynolds ceilings and sampled
frequency bands remain attached. A cold native sparse LDLT prepares the
steady response and checks the original momentum residual. Whole-interval
Reynolds admission uses passive energy contraction, including the steady
response residual. Accepted branch history can be imported/extracted and
continued under new viscosity without resetting its profile. Prepared
advancement reuses its state, flow and heat buffers; cancellation or refusal
preserves accepted history and reports.

The complete state/port augmentation is bounded by
`MAX_LAMINAR_NETWORK_AUGMENTED`; no original memory arm is truncated to fit.
Each connected component needs a pressure node, which may serve only as
the gauge of a closed loop. Junctions have zero storage and no modeled
entrance, mixing, area-change, thermal, cavitation or compressibility loss.
A supplied pressure head is an explicit physical source, not an inferred
pump constitutive law. Existing native pipe model and band limitations apply.

G1/G3 tests cover original single-pipe/series limits, independent
Poiseuille/Kirchhoff mesh equilibrium, loop head work, separate gauges,
viscosity/history continuation and refusal atomicity. The source-level
continuity tests cover topology, kinetic orthonormality, closed cycles,
dead ends, disconnected gauges and unequal inertias.

### Prescribed delivered flow through original network ports

`LaminarNetwork::prepare_flow_control` solves held pressure increments at selected
original reservoirs or directed branch-head sources for requested interval-mean
flows. Delivered volume equals that mean times the prepared duration; the
instantaneous flow remains the full original inertial/memory trajectory. The
coupled mobility comes from the cached native mean-output map, including every
initial profile coordinate. Cold preparation checks independence and scaled
positive response and factors once with native Cholesky. Warm calls reuse all
buffers, perform no new factorization or trajectory probes, and refine against
the complete original response rather than a fitted or symmetrized replacement.

The original component gauges are retained even when a reference reservoir is
controlled. Accepted pressure increments and gauged port pressures are reported
separately so a small physical drive is not rounded away beside a huge common
absolute pressure. Actual integrated port flows must satisfy the caller's
absolute/relative tolerance. The same complete trajectory, per-branch heat,
source work, continuity and whole-interval Reynolds checks gate acceptance;
state, clock, required pressures, flows and reports publish together. Dependent
controls, unresolved response, unsupported regime and cancellation refuse
without changing accepted results. This is a held-pressure volume-delivery
inverse of the existing incompressible developed-Newtonian network, with its
original rigid-wall, frequency-band and dense-size limits; it does not impose
instantaneous flow or introduce pump/compliance/cavitation physics.

## Unsteady laminar branches between common reservoirs

`graph::inertial::laminar::parallel::LaminarParallel` connects original
`LaminarPipe` owners between the same two held pressure reservoirs. Each branch
retains its actual material, geometry, viscous-profile history and Reynolds
ceiling. Total flow is the sum of the original trajectories at every instant;
the flow split is not frozen to steady Poiseuille conductances. Different
branches can carry opposing flows during a reversal.

Prepared advancement reuses each original exact propagator and independent
heat integral. Whole-bundle reports add actual pressure work, mean/profile
kinetic storage and irreversible heat; `branch_steps` retains each original
numerical admission and local heat. No joint solver certificate is invented.
Construction and preparation allocate; warm advancement reuses branch buffers
and publishes all histories and reports together only after complete admission.
`with_sources` preserves the original profiles and common clock under each
branch's admitted viscosity continuation, with density and geometry unchanged.

The summed condensed relation can solve a held pressure for an authored
**interval-average** total flow. It does not prescribe instantaneous flow or
model a sealed junction. Frequency response sums original admittances and
retains every branch's original band. The common density, developed Newtonian
flow, rigid no-slip walls and sampled-frequency approximation of `LaminarPipe`
remain; reservoir storage, mixing, manifolds and entrance losses are not inferred.
G1/G3/G4 tests compare unequal-radius startup/reversal to the independent radial
Stokes series, symmetric branches to the original pipe, zero-drive dissipation of
kinetic storage, signed total-flow pressure control, viscosity continuation,
harmonic parallel addition and atomic failure of a later branch.

## Geometry-derived laminar duct resistance

`MaterialFlowResistance::circular_laminar` and `rectangular_laminar` compile
straight, stationary, no-slip Newtonian duct dimensions into the existing
material-dependent pressure-loss law. Circular resistance uses Hagen-Poiseuille;
rectangular resistance retains the full Fourier side-wall correction, including
the square-duct and parallel-plate limits. Hydraulic diameter qualifies Reynolds
number but never substitutes circular drag for a rectangular cross section.

Each material evaluation supplies its actual dynamic viscosity and density.
Existing signed flow, pressure inversion, graph/thermal coupling, tangents and
dissipated-power accounting retain their owners. Six deterministic exponential
evaluations compile rectangular geometry once; no series or allocation enters
subsequent pressure/heat trials. The accelerated series' omitted correction is
below 1e-23 before floating-point roundoff, not a physical error certificate.

`annular_laminar` admits concentric stationary no-slip cylinders and compiles
the integrated Newtonian radial Stokes solution into the same owner. Annular
area supplies velocity and twice the radial gap supplies Reynolds length;
circular drag is never inferred from that diameter. A positive 17-term series
evaluated only at construction retains the cubic thin-gap conductance without
cancelling fourth powers. A vanishing but nonzero core retains its logarithmic
drag; exactly zero inner radius selects the original circular law. Eccentricity,
moving walls and annular thermal boundary conditions require separate models.
The pressure/convection regression places the annulus in series with a circular
heated section, so its changed viscosity changes both flow and downstream heat
transfer without applying a circular heat correlation to annular geometry.

The caller selects a Reynolds ceiling in `(0,2300]`; this is an applicability
screen, not a transition prediction. A developed Newtonian profile remains an
authored assumption. Entrance effects, non-Newtonian stress, moving/slipping
walls and unsteady-profile corrections are outside these constructors. Gas use
retains the existing frozen-density approximation. Focused regressions compare
analytic circular flow and an independent rectangular series, side exchange,
SI rescaling, signed dissipation and the original material/coolant heat path.

## Reynolds-dependent turbulent pipe resistance

`MaterialFlowResistance::circular_turbulent` explicitly selects developed
Newtonian circular-pipe Darcy-Weisbach loss with the implicit Colebrook-White
factor. Length, diameter, equivalent sand-grain roughness and a Reynolds band
are authored inputs. The band must lie within `[4000, 1e8]` and roughness/D
within `[0, 0.05]`; transition, entrance, non-Newtonian and compressible effects
are not inferred. Roughness is not substituted from a material label or RMS
surface height. The formula follows NASA KSC-STD-Z-0017 Rev A, §7.4.1, Eq. 4
and §7.5 (Colebrook 1939, DOI 10.1680/ijoti.1939.13150).

Each trial reuses the original current density and viscosity for Re, signed
pressure and `dp*Q`. A bounded allocation-free scalar Newton solve computes
the friction factor and its implicit derivative together; the reported
mass-flow tangent includes the change of friction with Reynolds number.
Geometry/calibration is retained in the evaluation rather than replaced by a
frozen quadratic fit. This owner feeds the existing liquid pressure, parallel,
pump and optional viscous-heating paths. Hydraulic heat enters the existing
enthalpy trajectory exactly once only when that dissipative path is selected.
The new pressure law does not supply a rough-wall heat-transfer correlation.
Thermal correlations remain separate declarations; a smooth-tube thermal
correlation still requires a smooth tube and its own supported flow regime.

Frozen-state `SeriesPressureFlow` keeps its analytic inverse for constant
linear/quadratic laws. Paths containing this turbulent law instead solve
inside the intersection of every original section's Reynolds support, with
one material snapshot per source, exact pressure tangents, cancellation polls
and at most 96 scalar iterations. Reverse flow uses the same passive curve;
no supported root means refusal. Fixed-coefficient graph lowering rejects
this mode rather than freezing its friction factor. This remains a
quasi-steady engineering closure, not a resolved turbulence calculation.

## Static and total fan shaft power

`qoi::FanPowerSpec` requires an explicit `FanPowerBasis`. Static efficiency
uses native AMCA static pressure times flow. Total efficiency adds the outlet
kinetic power of every original fan before dividing by efficiency; series
members carry the complete flow and parallel members carry flow/count.
The caller supplies actual density, physical outlet area per fan, count and
arrangement. Monotone positive-flow envelope propagation includes the cubic
kinetic term, and basis/geometry are bound into the existing QoI identity.
The same result flows through registered extraction and candidate QoI rows.

The result is mechanical shaft input. Motor/drive electrical loss, heat
deposition, geometry/density uncertainty, and total-pressure operating-point
rematching are not inferred. Supplied efficiency must apply to the declared
operating range; this is not a no-load motor or shutoff-power model. The native
static-pressure fan/network approximation and estimated evidence remain.

## Material-dependent enclosure heat in stored cooling

`StoredContactRadiation::MaterialClosedEnclosure` and
`MaterialOpenEnclosure` compose original source-dependent gray emissivity
with the stored solid, temperature-dependent contact, pump and finite coolant
transaction. Every mapped boundary receives its exact prepared trace owner.
Each joint h endpoint refreshes nodal emissivity and the actual reflected
irradiation system through the original conduction Newton solve.

The cooling adapter preserves source/card/chart identity and work allowances.
It returns only coolant convection to the fluid stream. Closed surface
exchange remains internal; the signed physical-surface sum from an explicit
black bath enters the solid external account once. Rejected thermal,
hydraulic, source or cancellation trials publish neither history.

The G1/G3 coupled regression checks closed and open cases against exact P1
moments of epsilon(T)*T^4 and independent two-surface radiosity elimination,
then checks the combined solid/tank first law and coolant return energy.
The current static geometry, patch-uniform irradiation, quasi-steady flow and
single-phase coolant assumptions remain. Native integration execution is
pending the shared thermal dependency build.


## Hydraulic and wall feedback on transient convection

With `transient-transport`, `transport::transient::coupled::advance_quasi_steady_graph`
resolves the native hydraulic graph and optional wall coefficients at every actual
thermal endpoint trial. Native Newton-Krylov and fallible directional differences
solve the original backward-Euler cell-energy residual. Hydraulic and heat payloads
are returned from an independent final evaluation at the accepted enthalpy field;
source refusals, cancellation and bounded endpoint-evaluation work remain active.

The existing thermal owner retains cell mass, caloric inversion, latent storage,
event boundaries, authored sources and independent local/global joule gates.
Explicit external mass streams supply their native enthalpy sources once; internal
upwind transfers remain conservative. The callback supplies only nonadvective
source/link/boundary coefficients and must be repeatable during derivative probes.
Hydraulic topology and branch identities are fixed. Original source applicability,
reference-density storage admission and complete cellwise mass balance still apply.

This is quasi-steady momentum with finite thermal inventory. Reversals update the
endpoint upwind direction; locating a within-step reversal, fluid inertance and
compressible storage require their respective owners. Omitted dissipative and
anomalous buoyancy work have explicit watt and joule allowances. Algebraic closure
does not replace time/control-volume refinement or experimental validation.

## Material-dependent graph drag and wall exchange

The natural-circulation material-response adapter rebuilds native graph loss
coefficients from explicit MaterialFlowResistance geometry and original material
states. Linear drag uses source viscosity; quadratic inertia and conserved flow
use the declared Boussinesq reference density. The original resistance owner
still admits Reynolds number, tangent and dissipated power. Each exchange uses
its explicit geometry and the original forced-convection correlation, retaining
material owners, Nusselt evidence and h. Source or correlation refusal is terminal;
constant-heat-flux cards cannot silently drive a constant-wall-temperature channel.
Topology, branch identities and segment ordering are explicit. Bulk/film sampling,
geometry/calibration support and material uncertainty remain caller-owned. This
binding is not a new constitutive database or a universal pipe-friction model.

## Coupled material response in steady natural circulation

natural_circulation::response::solve_circulation_with_response closes changes in
native drag and wall coefficients around the existing density/gravity/flow/heat
solve. The response callback receives actual branch means, declared-order segment
means and signed flows, and retains its material/correlation owners in the result.
Topology, exchange regions and wetted areas cannot change. Acceptance independently
checks unrelaxed temperature and flow residuals plus fresh linear/quadratic drag
and hA coefficients evaluated at the solved field. The inner buoyancy solve retains
its separate unrelaxed density and hydraulic-confirmation gates.

The previous integrated density profile warm-starts the next inner solve. Source
response calls, buoyancy samples and hydraulic solves have explicit total bounds.
Optional IQN-ILS works in caller-scaled coordinates; a direction change clears its
history and adopts the solved state without interpolating through artificial zero
flow. Reference density and Cp remain fixed conservation quantities. Spatial
property resolution and correlation bulk/film conventions remain explicit; no
transient inertia, coupled adjoint, stability or experimental-validation claim.

## Simultaneous contact and radiation in stored solid cooling

`stored_solid::contact::newton` adds
`advance_stored_solid_with_pump_and_contact_newton` to the original
immutable-history solid/tank/pump transaction. Its explicit
`StoredContactRadiation` selection admits no radiation, pointwise ambient
exchange with fixed or original-material emissivity, a closed gray enclosure,
or an open enclosure with one black reservoir. Each actual coolant-boundary mapping receives its original pointwise
owner; the complete contact, conductivity, storage and radiation tangent is
solved together by `TemperatureContacts`. Original heterogeneous enthalpy
charts, latent storage and source-domain refusals remain active.

`StoredContactRadiation::MaterialAmbient` retains the original material patch
bindings and exact storage caloric objects. Each actual coolant trial prepares
the material-radiation owner on its mapped boundary; current nodal temperature
and liquid fraction refresh the pinned emissivity within the same contact Newton
solve. Original convection alone heats the coolant return, while signed radiative
bath heat enters the external account once. A focused regression distinguishes
the current material response from its preparation snapshot using independent
P1 fourth/fifth moments and the complete solid-plus-tank first law.

`ContactNewtonStoredCoolingControls.thermal.solid` is the sole step policy
and must satisfy the contact owner's Armijo, zero-relative-tolerance and
absolute joule-target admission. One joint enthalpy solve is counted per
thermal callback against a whole-step cap; no nested contact/radiation
iteration product is introduced. The existing hydraulic, return-temperature
and wall-feedback owners still close the surrounding cooling circuit.

Coolant return energy uses only fresh original convective fluxes on the
network traces. Internal contact heat and reciprocal closed-enclosure exchange
are not external inputs. The signed physical-surface radiative sum, including
black-bath absorption when negative, enters the external heat account once.
Existing hydraulic heat policy, regional watt gates, independent solid/tank
joule checks and atomic publication remain unchanged.

Authored regressions compare surface emission to exact P1 fourth moments
and two-surface reflection to an independent 2-by-2 radiosity calculation,
alongside unequal-temperature contact/coolant storage, constant-resistance
agreement and transactional refusals. Geometry and quasi-steady pipe transport
remain fixed. Ambient emissivity uses the explicitly chosen snapshot or original
material response; enclosure emissivities remain snapshots. These methods add no
fluid phase-change, moving-contact or empirical-validation claim.

## Transient graph convection with native thermal storage

With `transient-transport`, `graph::thermal::transport::transient` binds signed
hydraulic flows to explicit native enthalpy cells through `GraphEnthalpyTransport`.
Those cells retain their mass, volume, caloric law, latent heat and heat sources.
Reference density converts volume flow to mass flow; every storage node must
balance mass, including pressure datums. Inlet/outlet enthalpy streams are
explicit. The original thermal stepper owns combined conduction/advection CFL,
material-domain checks and energy accounting. Closed adiabatic circulation
needs no artificial reservoir because initial stored enthalpy fixes its state.
A bounded backward-Euler specialization solves advection plus fixed watts on
at most 64 cells, including the native enthalpy-datum correction, and accepts
only through the original local/global energy residual. General wall/contact
exchange remains available through the native explicit observer and implicit
endpoint-coefficient interfaces. Flow and density are frozen; cells are perfectly
mixed, and pipe residence-time resolution requires explicit spatial inventory.

## Native graph temperature–buoyancy feedback

`graph::thermal::natural_circulation::solve_natural_circulation` closes steady
temperature -> material density -> Boussinesq head -> hydraulic flow -> heat
transport on closed native graphs. Straight branches carry explicit axial
segment fractions; density is integrated along the actual exponential fluid
temperature profile, with ordering reversed with flow. Material validity and
density-departure bounds apply at both endpoints and interior samples. Paired
quadrature, source-query and nonlinear-iteration budgets are explicit.
Acceptance requires the unrelaxed density residual and a fresh hydraulic solve
at the resulting thermal densities to satisfy the flow-change tolerance.
Optional bounded IQN-ILS is acceleration only and resets after flow reversal.
The result retains heat, density, gravity-head, flow and omitted-dissipation
diagnostics. Pressure datums are not thermal supplies. Reference density/Cp,
drag and hA remain fixed; stagnant branches, unsupported material states and
excess omitted viscous heating refuse. No startup, compressibility, nonlinear
stability, coupled-adjoint or experimental-validation claim is made.

## Finite source-backed gas circulation

`feedback::transient::FiniteGasLoop` uses the same private finite-inventory
momentum/thermal transaction as `FiniteLiquidLoop`; the liquid API and native
physical owners are unchanged. Its exact registered `GasEnthalpyLaw` owns
NASA-9 enthalpy, bounded inversion, caloric stability and formation-energy
references. Native cells store M*h and signed upwind links carry original donor
h. Neither nonlinear heat capacity nor nonzero formation enthalpy is replaced
with Cp(T)*T or temperature averaging.

Reference density is evaluated at the caloric law's original context. Actual
gas snapshots reevaluate only temperature, preserving pressure, composition,
transport selections, source-pack identities and other source conditions.
Mass, volume, continuity and inertance remain fixed reference quantities;
source density departures act only in gravity. The original common density,
Reynolds, omitted-work, energy, momentum and rollback gates still apply.
An explicit Mach limit in (0,0.3] also admits the largest old/new/average speed
using each endpoint source's sound speed. A small average during reversal
cannot bypass either Mach or Reynolds admission.

This lane describes a frozen-composition, fixed-pressure, low-Mach Boussinesq
gas circuit. It does not evolve pressure, compressible mass, chemistry or
condensation, or infer a zero-flow wall film. Original source/caloric support
and all bounds are mandatory; an unsupported gas is never replaced by air.
Returned hydraulic source snapshots retain the original thermochemical
receipts. The native thermal ledger measures the original stored enthalpy and
actual heat sources independently of the separate bounded mechanical account.

`tests/finite_gas_loop.rs` uses a separately reduced nonlinear NASA polynomial
and cyclic enthalpy balance to check finite storage, donor transport and buoyant
acceleration, plus a reversal whose endpoint Mach limit must refuse atomically
although its signed-average speed is smaller. These are regression definitions,
not a claim of executed native tests or experimental validation.


## Finite pipe momentum and thermal inventory

`LiquidHeatPath` remains a quasi-steady, zero-pipe-storage march. The separate
`feedback::transient::FiniteLiquidLoop` binds actual native thermal cells and
their original liquid caloric/transport claims to the existing PHS inertance.
It admits a closed connected path, fixed reference density rho0, V=A*length and
M=rho0*V. Temperature-dependent source density supplies only gravity; source
viscosity supplies the existing loss law, with rho0 used for continuity, loss
inertia and pipe momentum. Every trial enforces the declared density-departure
limit and original source/Reynolds support, including Reynolds zero.

The original PHS discrete-gradient step advances momentum with
L=sum(rho0*length/area), H=p²/(2L). The original thermal backward-Euler residual
and signed upwind `Advected` transfers evolve the same finite pipe enthalpy.
Newton--Krylov closes their common endpoint; it does not own another physical
timeline. Heat transport uses the actual PHS step-average volume flow, including
zero and negative values. Its sign may differ from endpoint velocity during a
reversing step. Original thermal sources, links and boundaries remain active at
rest; no singular positive-flow correlation or invented zero-flow film is used.

Acceptance gates original-law impulse, kinetic-work, every-cell thermal and
total thermal residuals. The PHS coefficient trial is checked again against the
original nonlinear losses at its actual resulting flow. Source, callback, work,
solver or cancellation refusal publishes neither thermal state nor momentum.
Exact original caloric-object binding, geometry, mass and accepted time are
retained. Pure heat callbacks may be evaluated repeatedly; only their accepted
payload is returned. A shared budget counts all residual and derivative trials.

This is an explicit low-Eckert Boussinesq model, not a compressible total-energy
closure. Mechanical and thermal accounts stay separate. Hydraulic loss and
absolute local density-anomaly gravity exchange are bounded together in watts
and step joules using max(|Q_old|, |Q_new|, |Q_average|), without cancellation or
invented pressure-work heating. The original monotone loss law is evaluated at
that envelope speed with the same endpoint-held source states, also admitting
its Reynolds support. A symmetric reversal cannot hide omitted work in zero
average flow. These bounds describe the admitted discrete interval/coefficient
approximation, not unresolved continuous-time extrema or temporal accuracy.
Envelope losses stay separate from native signed-average discrete work. The
native kinetic receipt reports pressure supply, buoyancy work and dissipation;
the thermal receipt reports only its actual authored heat. A dt*R/L limit at
most two prevents an unforced midpoint-drag sign flip. Algebraic acceptance
does not certify time or spatial accuracy, crossing-event timing, boiling,
turbulence, evolving volume or a wall-film correlation.

Tests in `tests/liquid_pressure/transient.rs` use independent cyclic
finite-volume formulas and scalar inertance balances for heating from rest,
finite pipe thermal memory, step-average versus endpoint direction, signed
reversal, original viscosity response, a symmetric-reversal omitted-work cap,
and atomic failure/approximation limits.


`MaterialFlowPath::march_liquid_enthalpy_with_dissipation` integrates the original
resistance owner's `dp*Q` as uniformly distributed fluid heat in each region,
using the existing piecewise caloric exponential march. It requires constant
source density and equal flow areas; choosing it declares equal elevation,
unchanged bulk kinetic energy and complete local conversion of the specified
pressure loss. The source Cp chart remains authoritative. Actual wall heat is
`mdot*(h_out-h_in) - P_loss`, and only that wall heat enters the original Robin
row. The report separately retains hydraulic evaluations, hydraulic heat and
actual endpoint advected enthalpy gain; a compensated balance residual exposes
reporting error instead of defining the enthalpy gain by closure. It does not
model gas throttling, pump inefficiency or a
porous skeleton's share of dissipation. Unsupported source/Reynolds/geometry
states and cancellation return no partial trajectory. Analytic constant-Cp,
variable-Cp subdivision and original-owner refusal tests cover this opt-in path.
The shared caloric marcher refuses a nonzero physical exponential transfer if
division by mass flow underflows or its specific-enthalpy addition rounds back
to the old value. This applies to ordinary and dissipative marches alike:
unrepresented storage must not masquerade as wall cooling. A source-backed
large-reference-enthalpy regression covers hydraulic, wall-heating and
wall-cooling refusal while the same increments remain valid at an ordinary datum.
The hydraulic source shift also refuses a nonzero driving-temperature increment
that underflows or disappears on addition, before it can bypass that caloric guard.
The shared marcher also checks the accepted caloric temperature against the
original exponential's inlet/driving-temperature interval, allowing temperature
interpolation roundoff only. A nonzero enthalpy increment that rounds past the
wall therefore refuses instead of publishing an impossible heating or cooling
endpoint with a self-consistent but physically wrong Robin heat report.

`LiquidHeatPath::march_with_dissipation` closes regional outlet temperature,
convection coefficient and hydraulic power together. The original liquid law
evaluates viscosity at the arithmetic bulk mean and convection at the selected
bulk/film convention. Bounded outlet predictions establish a local bracket,
then bisection closes it; the solver never probes an unrelated caloric maximum
outside a narrower transport source. Every trial still obeys original source
and Reynolds validity. The final re-evaluation checks temperature, HTC and
hydraulic power with the caller's existing tolerances. Wall-only Robin rows
retain the actual wall load and source-backed endpoint states independently
report advected energy. `solve_pressure_driven_with_dissipation` composes this
trajectory with the existing bounded pressure root and requires zero omitted
hydraulic-power allowance. Analytic affine-viscosity, pressure-selected flow,
narrow-source and minimum-iteration regressions exercise the complete path.
`solve_pressure_coupled_iqn_with_dissipation` feeds that same heated trajectory
through the original solid IQN wall iteration and pressure root. Every sent
Robin row and every re-evaluation at the actual returned solid walls includes
hydraulic deposition, while regional solid balance compares wall heat alone.
The pressure root retains the hydraulic state from the accepted re-evaluation,
with one shared callback limit across all pressure trials. Source-density,
geometry, cancellation, regional power and zero-omission gates remain intact.
`tests/liquid_pressure/dissipative_solid.rs` covers the native conduction FEM
with separate source/hydraulic heat, temperature-dependent bulk/film properties,
and global-work, local-balance, density and omission refusals.
`ParallelLiquidNetwork::with_hydraulic_dissipation` selects this same branch
marcher once for pressure redistribution, common-solid IQN, pump matching and
finite-reservoir continuation. Original material/chart identities and the
selection survive every inlet rebind. Mixing and reservoir return use actual
mass-weighted endpoint enthalpy changes; wall heat and deposited hydraulic power
remain separate. The native stored-solid transaction counts that hydraulic heat
once as a combined-system input, while its internal solid/coolant gate compares
wall exchange alone. Its existing contact/radiation variants share this account.
Branch, network and reservoir omission allowances must all be zero in the
selected model; original local, global-work, source and atomic-publication gates
still apply. No extra pump-efficiency heat is inferred. Analytic parallel mixing,
variable-material feedback and two-step native solid/reservoir regressions live
in `tests/liquid_parallel/hydraulic_heat.rs` and
`tests/stored_cooling/hydraulic_heat.rs`. Ordinary networks retain their explicit
omission budgets and original omitted-heating behavior.

> Status: ACTIVE for beads `frankensim-extreal-program-f85xj.5.5` and
> `frankensim-extreal-program-f85xj.5.7` (conjugate coupling).

## Purpose and layer

`fs-airflow` is the L3 low-cost enclosure-airflow rung and the dependency-safe
consumer seam for the thermal vertical's decision-facing QoIs. It retains typed fan
pressure/volume-flow data, validates monotone interpolation, applies bounded
fan-law speed scaling, composes quadratic loss elements in series and parallel,
requires an explicit leakage branch, solves the fan/system operating point, and
combines that record with a steady `fs-conduction` solution without making either
lower-layer solver depend on the other.

The `conjugate` module closes that combination into a two-way exchange. It owns
the stream-wise air state the flow-network types do not carry, and the
partitioned fixed point between a solid conduction solve and that air path.
`fs-couple` has no coupling-loop driver — its only iteration entry points are
hardwired to a private added-mass fixture map, and its contract assigns the
driver to the consumer — so this crate is that consumer. Interface quasi-Newton
(IQN-ILS) does not exist in `fs-couple` and is not called.

Runtime dependencies are `fs-blake3`, `fs-conduction`, `fs-convection`,
`fs-couple`, `fs-evidence`, `fs-exec`, `fs-ivl`, `fs-ladder`, `fs-math`,
`fs-qty`, and `fs-regime`. Results flow outward as evidence-bearing typed
quantities, as a Re/Pr handoff to the existing convection rung, and as the
CHT ladder's correlation-rung transfer.

### Material-driven natural and forced finite-cell convection

`medium::convection::deposition::evaluate_material_cell_convection` binds either
an existing forced-convection card or an explicitly oriented natural-convection
geometry (vertical plate, horizontal cylinder, horizontal plate) to
the original finite-cell enthalpy owner. These laws share one source, exposed-area
and contact-coverage accounting path. Natural convection reevaluates Ra and h
at the proposed wall temperature inside the thermal residual; density, Cp,
viscosity and conductivity come from one admitted material transport snapshot.
Reservoir properties and explicitly supplied volumetric expansion remain frozen
for the interval. No liquid expansion or film state is inferred.

Exact thermal equilibrium returns zero power while retaining the material
snapshot, with no fabricated Ra=0 correlation result. Every nonzero contrast
retains the original card's domain restrictions. Repeated bindings, authored
heat-transfer overlap, invalid sources, work limits and cancellation refuse
transactionally. The existing forced-only entry point delegates to the same
ownership loop. `tests/material_cell_convection.rs` covers independent natural
heat-transfer calculations, material substitution, both heat directions,
negative expansion, nonlinear accepted endpoints, ambient heater startup,
mixed-law ownership and thermal-energy balance.

### Source-resolved film convection on the same finite-cell residual

The shared `AidingVerticalPlate` law supplies Manning-Qureshi average laminar
mixed convection on an unconfined isothermal vertical plate. Its source-backed
cell/link and `FilmGasHeatLaw` adapters use one mean-film property query and the
full wall/bulk contrast. Film averaging is a declared frozen-property adapter
approximation; it is not attributed to the paper's constant-property theory.
Signed vertical speed must aid the actual signed expansion response. The caller
also supplies an explicit maximum `|beta delta(T)|` below one. Opposing flow,
excessive density departure and unsupported laminar/Prandtl groups refuse; no
cube-difference, invented conduction floor or turbulent extension is supplied.
The same card reaches its natural limit continuously as fan speed tends to zero.
Plate height, masked exchange area, original source receipts and native thermal
owners remain separate. No mass/energy store or thermal solver is added.
Focused `mixed_convection` and finite-gas film tests cover actual signed heat,
source state/context, stopped/running flow and reciprocal energy accounting.

`medium::convection::deposition::source::evaluate_source_cell_convection`
re-evaluates the original fluid owner at an explicitly selected reservoir or
mean film temperature for every proposed thermal endpoint. The property
temperature never replaces the actual reservoir temperature in buoyancy or
signed heat flow. Density, heat capacity, conductivity, viscosity and known
signed isobaric expansion share one checked source state. Natural convection
requires sourced expansion; forced convection does not invent it when absent.
The `Natural` horizontal-plate, cylinder and sphere laws require `FilmMean` and
refuse bulk-property substitution before querying sources. Signed expansion and the
full wall/bulk contrast select plate stability; exposed contact-masked area
does not replace the explicitly authored full-plate `A/P` length. The same law
resolution applies to reciprocal finite links, so the native thermal owner
debits and credits one heat transfer. Focused `natural_geometry` tests check
orientation, negative expansion, retained pressure/composition, independent
film-property heat calculations and closed-pair energy conservation.

The result retains each concrete material owner and its thermal snapshot beside
the original correlation, exposed-area and heat-source results. Prepared source
laws retain pressure, composition and selected claims across trials; the generic
transport interface does not establish that identity on a caller's behalf.
Source support, duplicate-boundary ownership, work limits and cancellation
remain transactional. The prescribed reservoir has no simulated finite energy
store, and each correlation still requires its own supported property convention.
`tests/material_cell_convection/film_source.rs` checks independent film-property
algebra, signed expansion, source-only viscosity substitutions in accepted
backward-Euler finite-body balances, retained claims and atomic refusals.

### Original liquid convection on implicit conductive surfaces

`medium::convection::implicit::ImplicitLiquidConvection` now implements the
native `ImplicitGrayExteriorLaw` on the same bulk caloric object consumed by
`IntegratedCaloricGrayBoundary::with_exterior_law`. The existing owner closes
the original k(T) half-cell, optional coating, convection and gray radiation
at one actual exterior temperature. It retains separate bath powers, energy,
entropy and the original implicit heat timeline; this adapter adds no thermal
store, material resistance or second surface solve.

Preparation retains the original liquid card, selected claims, fixed pressure/
composition and source receipts. Supported laws are original linear temperature
curves and constant scalars with explicit thermal validity. Each native segment
is prepared once; its derivative is taken from original source coefficients,
not differences between clipped endpoint queries. Surface trials use bounded
binary searches and the existing evidence-bearing correlation evaluator,
without re-querying the material database. Sampled source surfaces, discrete
sources, augmented bundles and unsupported spans refuse rather than being
fitted into replacement temperature laws.

Mean-film transport includes the factor one half in every temperature chain.
Sieder-Tate and Whitaker retain reservoir Re/Pr/k and differentiate only their
same-source wall viscosity. A kinematic-viscosity source keeps the complete
rho(T)*nu(T) derivative. Natural geometries retain signed original expansion
(from the admitted density path or selected beta claim), full surface/bath
contrast and the original stable/unstable selector. Moving-input partials,
source kinks and empirical regime endpoints are not guessed. Dittus-Boelter
also checks its declared heating/cooling convention against the actual drive.
Laminar aiding mixed convection uses the original Manning–Qureshi card and
analytic Re/Ra/Pr partials on the same prepared film source. Signed expansion,
signed stream speed and the explicit density-departure bound pass through the
shared material admission; opposing flow and unsupported source/regime states
refuse. Pure natural flow keeps Re identically zero and requires no Reynolds
derivative. Finite plate strips are not exposed by this adapter.

The native flux-state override admits the natural sphere's original Ra=0,
Nu=2 equal-temperature limit with total heat tangent 2k/D, without claiming a
finite derivative of the singular coefficient. Identically zero driving from
zero gravity or locally constant zero expansion retains the same sphere
conduction branch. Other natural geometries retain their positive card minima.
At exact mixed-flow thermal equilibrium the native heat tangent is the original
pure-forced coefficient, even when the buoyant coefficient derivative is
singular. At exact resting equilibrium both heat and its tangent vanish as the
original laminar natural component's limiting heat response; no Nusselt evidence
or conductive floor is invented. Original Pr/source admission still applies.
Fluid and exterior support intervals remain distinct; actual bulk/film and wall-viscosity
property coordinates are admitted on each trial. The targeted
`material_cell_convection::film_source::implicit` regressions cover analytic
film and wall-viscosity chains, signed natural/mixed geometry, pure-flow limits,
source refusals, sphere/mixed equilibrium, integrated half-cell heat and native backward-Euler
energy and separate-bath entropy accounts.

### Original gas convection on the same implicit conductive surface

`medium::convection::implicit::ImplicitGasConvection` now binds the original
`IsobaricGasLaw` to that native `IntegratedCaloricGrayBoundary` interface. It
retains the source pressure, fixed composition, original transport choices and
phase admission while resolving convection at the actual exterior temperature
after the body's k(T), half-cell and coating resistance. Prescribed gas and
radiative baths remain separate external heat/entropy destinations; this
adapter creates neither a finite gas store nor another thermal integrator.

Gas and liquid owners share one correlation/chain-rule evaluator, geometry
admission and equilibrium heat-tangent implementation. The common public
`ExteriorConvectionSettings`, `ExteriorConvectionResponse` and
`ExteriorConvectionError` retain the original `LiquidExteriorConvection*`
names as compatible re-exports. Mean-film laws include the one-half temperature
chain in rho, Cp, k and viscosity. Sieder-Tate and Whitaker query the SAME gas
at the wall solely for the original viscosity correction, retaining bulk
Re/Pr/k. Natural and aiding mixed convection retain the original isobaric
ideal-gas expansion 1/T and its derivative -1/T^2, actual signed drive,
source-domain checks and exact sphere/mixed equilibrium heat limits.

Analytic source tangents come from the original native Sutherland, USSA/Eucken,
named Wilke/WMS and supported thermochemical gas models, including actual
NASA heat-capacity slopes and composition-dependent mixing denominator
derivatives. The thermochemical path differentiates its selected original
transport claims, including original temperature/pressure surfaces at fixed
pressure. Active film/wall adapters retain prepared source geometry once;
pressure edges and inward thermal endpoints keep their supported derivatives.
Actual trial queries preserve original source selection and validity, while
unsupported derivative forms, holes and genuine source kinks refuse.
No source fits or production finite differences are introduced. Bulk-held
circular developing laws retain exactly zero wall-property derivatives and do
not require an active transport tangent. Each bulk,
film and viscosity query rechecks declared temperature support. Original
humid-gas saturation and named-species domains remain mandatory; rectangular
four-isothermal-wall sections still require a whole-perimeter thermal owner
and cannot enter a pointwise exterior boundary. Analytic derivatives do not
supply an empirical model-accuracy bound or expand any correlation's domain.

Focused gas regressions exercise original pressure/composition, exact material
queries, film and wall-viscosity chain rules, natural/mixed expansion, phase and
geometry refusal, and bulk-held and active thermochemical source support. A
manufactured NASA/transport homogeneity identity independently checks the
nonzero Cp derivative through the film-temperature correlation. A native
conducting-wall regression independently balances the analytic k(T) half-cell,
finite plate strip and gray radiation, then checks accepted backward-Euler
storage and separate bath energy/entropy accounts for dry air, a named gas
mixture and an original thermochemical mixture with variable heat capacity.
These runnable tests do not claim experimental validation.

`source::observer::SourceConvectionObserver` supplies these heat rates to the
original explicit `DepositionSimulation::advance_observed` timeline and composes
with other `DepositionObserver` owners. Every actual substep reads its starting
wall/film state and exposed area; accepted endpoint samples retain freshly
resolved material owners. Scheduled events invalidate the sample, and the next
substep uses the updated contact coverage and boundaries. The observer's signed
interval heat is a share of the existing external-heat audit, never a second
energy addition. Endpoint source refusal or cancellation discards both proposed
histories. Resolver cloning must preserve the native observer's rollback rules.

One observer may bind the complete progressively deposited assembly. Its total
binding count (active plus future) is bounded by `maximum_cells`; duplicate
IDs refuse before any source query. On the first source query of each advance,
unborn IDs must match native scheduled births. This schedule scan is checked
once, with cancellation points, and is not repeated on thermal substeps.
The lazy `new` constructor cannot validate unknown IDs during empty spans or
zero-time endpoints, where the native timeline performs no source query.
`new_for_timeline` admits the same complete roster immediately, with caller
checkpoints and no fluid-property queries, before either kind of advance.
Every new advance resets schedule admission against its actual thermal input;
the immutable binding index is reused, not a previous timeline's future roster.
Unborn cells contribute no heat, CFL restriction or fluid-property query.
The birth callback admits only that cell's actual material/correlation state,
including when an advance ends exactly at birth. Unsupported incoming states
discard the candidate advance. Same-time contacts may still change exposure:
the birth check neither applies provisional-area heat nor validates coverage
before those contacts arrive. Actual source and timestep evaluation follow all
same-time events and retain the original exposed-area and boundary-ownership
checks. Retained endpoint samples remain invalidated by events. Focused
`film_source::observer` tests cover staggered and simultaneous births, independent
piecewise heat/CFL balances, supported contact coverage, future-source domain
refusal at birth, roster budgets and cancellation.

The same start-state hA now tightens native subcycling through
`DepositionSimulation::limit_timestep_for_additional_conductances`, using the
original caloric slope bound and cell mass. Its additional decay rate composes
harmonically with the already admitted step; independent minima would reuse the
same monotonicity allowance. No heat is applied a second time. Pending
conductances are invalidated at events, acceptance and a new advance.
This adapter retains explicit staggering and time-refinement requirements;
the frozen-coefficient bound is not a general nonlinear-stability guarantee.
It does not evolve reservoir energy, change source policies or introduce
another thermal integrator. `tests/material_cell_convection/film_source/observer.rs` compares
independent finite-body updates with simultaneous native spatial conduction,
event-clipped exposure changes, retained endpoint sources, interval heat accounts,
late source-support refusal and cancellation.

### Material fluid cooling in the complete circuit/heat network

`source::heat::SourceConvectionHeatNetwork` resolves the existing prepared heat
network and source-backed fluid convection at the SAME trial enthalpy endpoint.
Its circuit entry calls the original contact-aware circuit/enthalpy driver.
Bulk and interface electrical losses, nonlinear conduction, source-integrated
layers, thermal contacts, enclosure radiation and fluid cooling therefore share
the original finite caloric cells. Changing a fluid property law can change
cooling, body temperature, material resistance and actual capacitor/current
evolution; no component-name dispatch or separate thermal step is introduced.

The fluid resolver retains its source models, fixed pressure/composition and
claim selections across trials. Reservoir or film property temperature comes
from the existing convection owner; exposed area and contact coverage stay
with the thermal owner. Duplicate authored ambient coefficients refuse.
All accepted bulk/interface/radiative and fluid source records remain available.
Convection is external heat from prescribed reservoirs; the closed enclosure's
signed net radiative defect is excluded from external supply by an additional
combined gate using the original total-energy tolerance.

This retains the native generally first-order endpoint coupling, source domains
and uniform-cell limits. Time refinement remains necessary. It does not evolve
reservoir inventory, flow, mass exchange, pressure work or moving geometry.
`tests/material_cell_convection/film_source/heat_network.rs` supplies an
independent temperature-dependent RC/film-cooling reference, fluid source
substitution, current/storage/energy checks and transactional refusals.

### Finite coolant in the complete material heat network

`SourceConvectionHeatNetwork::resolve_with_links` and
`step_candidate_with_links` compose finite-fluid convection with prescribed
reservoir convection and the same complete heat/circuit assembly. Initial
link conductances enter `MaterialHeatNetwork::prepare_with_link_conductances`;
each subsequent residual resolves them afresh from both actual endpoint
temperatures. The original cells retain their own caloric laws and fixed
inventories. Finite convection stays in reciprocal link conductances, never
in the external source-power account.

Fourier bindings and additional convection bindings must partition every
native link exactly. Missing, duplicate, inactive, nonfinite or competing
ownership refuses; ordinary preparation/resolution remains strict.
Integrated conduction may replace only its original Fourier links, and a
finite convective link cannot also carry a separate thermal-contact law.
The report retains the complete heat assembly, prescribed-reservoir sources
and finite-fluid source evaluations together.

The focused `heat_network/finite_links.rs` regression uses one liquid card
for its caloric storage and film transport. Viscosity and heat-capacity data
substitutions alter both finite temperatures and the original RC charge/current
solution, compared with an independent scalar balance. Internal exchange
cancels and capacitor loss supplies total thermal gain. The fixture also
checks source pins, link orientation, missing/conflicting owners and refused
source trials. Uniform cells, fixed geometry/mass, prescribed relative flow
and native first-order endpoint coupling remain the model limits; advection,
fluid momentum, pressure work and moving phase fronts need their own owners.

### Evolving gravity circulation in native finite coolant inventories

`source::heat::relaxation::pressure::gravity` connects the existing material
series-pressure solver to the original frame/circuit/enthalpy transaction.
`step_gravity_driven_relaxation_with_links` resolves actual bulk source density
and viscosity at every proposed endpoint, selects signed circulation from
imposed pressure plus gravity, derives the section velocities, and supplies
that same mass rate to the native conservative upwind enthalpy graph. Reciprocal
wall heat uses the original source convection and original finite-cell links.
The accepted payload retains the exact source states used by that endpoint;
there is no second coolant state, heat integrator or post-acceptance update.

`SeriesPressureFlow::solve_with_boussinesq_gravity` uses explicit
`BoussinesqDensity`: one reference density for continuity, volume conversion,
inertia and hydraulic Reynolds, with source density retained separately for
buoyancy. Each original resistance evaluation exposes both densities; viscosity
and thermal properties remain source-resolved. The loop requires common
reference mass/volume, distinct inventories and exactly closed connected
geometry. Every bulk density must satisfy the caller's relative-departure
bound. Uniform-density rest cancels exactly, and reversal uses the same signed
flow and native donor-enthalpy convention. Film properties remain an explicitly
separate convection evaluation, never the bulk hydrostatic density.

This is a quasi-steady low-Eckert model. Imposed-pressure power, buoyancy power,
local hydraulic loss and their independent mechanical residual are reported
separately. None is disguised as an external heater. The sum of their absolute
power magnitudes and its interval integral must satisfy caller-provided W and J
omission bounds, including during nonlinear trials; opposing drives cannot
hide omitted conversion through cancellation. These bounds cover the reported
hydraulic terms, not every error in the Boussinesq approximation. The original
thermal and circuit accounts close under that explicitly admitted omission.
This is not an exact total thermodynamic-energy, compressible pressure-work,
fluid-inertia or changing-volume solver.

Unsupported material/correlation states, density departures, work budgets and
either cancellation authority reject the complete mechanical/electrical/caloric
candidate. A stagnant forced-convection correlation still refuses if outside
its domain; no artificial minimum speed or heat-transfer floor is supplied.
Adiabatic loops admit exact rest. Focused G1/G3/G4 tests cover analytic reference
density hydraulics, temperature-driven circulation, reversal, unchanged stored
mass/geometry, native sensible-energy balance and whole-history rollback.

### Shared heat assemblies with material field storage

The `source::heat::assembly` module implements the fs-couple `HeatAssembly`
interface for prescribed fluid reservoirs and the borrowed `with_links` view
of finite-fluid convection. Both delegate to the existing complete source
resolver. The exact prepared heat owner, all fluid evaluations and original
heat-network sample remain accessible; no fluid type or material name selects
a physical driver.

`fs_couple::circuit_enthalpy::field_storage::heat::field_storage_candidate`
consumes that same assembly in the original C(T)/L(T), current and enthalpy
residual. Original field-parameter work, entropy and reversible heat remain
separate from Joule loss. Finite coolant exchanges reciprocal internal heat;
prescribed reservoirs supply external heat. Closed radiation's numerical net
is excluded by the original physical energy gate. Source errors remain typed
and native cancellation/work failures keep their original classifications.

The finite-links field-storage regression compares both capacitor and inductor
histories against independent midpoint and total-internal-energy equations.
Source-only viscosity, coolant Cp and field-temperature-slope substitutions
change accepted charge/flux and temperatures. It covers continued steps,
source identity, missing finite-link ownership and candidate refusal/retry.
Original whole-volume field capacity admission, uniform thermal cells, fixed
geometry and timestep-refinement requirements remain unchanged. Pressure,
composition and relative flow retain their original provider ownership.

The same assembly also feeds
`fs_couple::circuit_enthalpy::piezoelectric::heat::step_candidate`. Original
rod/circuit reciprocity, source-bound load resistance and finite coolant heat
share one accepted candidate. The base ceramic coefficients remain held;
cooling changes actual current and continuing motion through the original load.
Its finite-links regression compares independent midpoint and two-cell energy
equations under viscosity/Cp substitutions, preserves full source records,
and exercises late typed source refusal, circuit budget and cancellation.

The temperature-dependent ceramic entry
`fs_couple::circuit_enthalpy::piezoelectric::thermal::heat::step_candidate`
consumes this same provider and an explicit whole-ceramic thermal binding.
Compliance, signed coupling, permittivity, load resistance and fluid transport
resolve at one proposed endpoint. The original ceramic driver owns entropy,
internal energy and reversible heat; the native circuit owns Joule loss.
The complete coolant sample survives in `report.coupled.source.network`,
and typed provider refusals remain under `PiezoelectricThermalError::Heat`.
The finite-links ceramic regression checks independent constitutive, midpoint
and two-body internal-energy equations over continuing steps, source-only
coolant/ceramic substitutions and late source/budget/cancellation refusal.
Uniform fully solid material, fixed geometry and original capacity/source
limits remain; thermal strain, depoling and dielectric loss are not inferred.

Thermoelectric conversion and its existing reciprocal field-motion variant
also consume `HeatAssembly` through the free functions in
`fs_couple::circuit_enthalpy::contacts::heat::thermoelectric`. Finite coolant
and prescribed reservoirs enter the original Kelvin/circuit/enthalpy residual.
Integrated Seebeck voltage, terminal Peltier heat and interior Thomson
allocation retain their original material source and destinations; bulk and
contact Joule heat are still deposited once. Complete accepted fluid samples
remain in `report.coupled.source.network`, and provider errors remain typed
under `ThermoelectricNetworkError::Heat`. Closed radiation and finite coolant
exchange remain internal to the same energy account. The finite-links Kelvin
regression checks independent four-cell and charge equations, reversible
entropy, continued state and source-only viscosity/Cp substitutions. Existing
leg trace/contact restrictions, fixed heat geometry and first-order coupling
remain; the moving field/support materials keep their declared native scope.

### Material graph feedback with complete source-resolved heat

`fs_couple::constitutive_graph::enthalpy::heat_assembly_candidate` consumes
the same complete provider for public material-law feedback. Actual trial
temperatures drive both original graph inputs and coolant/material transport;
the existing graph and thermal residuals remain the sole state owners. Finite
coolant exchange is internal, graph watts are routed once, prescribed heat is
external and closed-radiation roundoff retains its separate account. The full
accepted fluid/heat sample remains in `report().source`. Prepared circuit
bindings and resolved resistor/contact work refuse so graph heat cannot omit
or duplicate another electrical owner. Typed source/cancellation failures
remain native; other provider refusals retain their type under `Transport`.

The finite-links graph regression checks a source-bound Ohmic law under a
solved voltage constraint against independent current and two-body heat
equations. Continued thermal states and viscosity, coolant Cp and conductor
resistivity substitutions use the same original solver; late source, budget,
cancellation and retry preserve accepted state. Original graph work/storage
obligations, source ranges, fixed geometry and first-order time accuracy
remain authoritative.

### Hysteretic core and winding feedback with finite coolant

`ThermalHystereticCircuit::step_candidate_with_heat_assembly` consumes this
same complete source-resolved heat provider. Its original current, magnetic
memory and finite-enthalpy solve resolves winding R(T) and coolant at each
trial endpoint. Winding Joule heat and core pinning heat enter their original
distinct cells once; reciprocal finite coolant exchange remains internal.

`bind_with_link_conductances` admits the initial finite-link conductances
through the original complete/disjoint material-link checks. It does not
cache those values as a cooling law: continued assembly trials resolve them
again, and a bare step still refuses links without their physical owner.
The exact bulk material network, full accepted fluid samples, source identity,
typed provider errors and original energy gates remain authoritative.

The focused finite-coolant fixture independently checks current, play memory,
core pinning heat and the combined energy balance over two intervals. Source
substitutions change coolant viscosity, coolant heat capacity and magnetic
pinning field. Missing link ownership, a later typed source refusal, work
exhaustion and cancellation preserve the accepted states. Magnetic fit
coefficients remain constant within their admitted thermal interval; the
integration adds no unmeasured thermal magnetic law.

### Finite coolant in reciprocal magnetostrictive actuation

`ThermalMagnetostrictiveCircuit::step_candidate_with_heat_assembly` uses the
same complete coolant provider for the original biased winding, rod motion
and finite thermal residual. Source-only coolant changes alter winding R(T),
actual current, linkage and motion. `new_with_link_conductances` admits the
initial finite links with the original complete/disjoint ownership checks;
each continued trial resolves those links again.

The original total-current Joule heat includes DC and cross terms, and the
storage account retains `I_bias * delta_lambda`. Neither term is replaced by
an incremental-only estimate. Finite exchange stays internal, accepted fluid
samples stay inspectable, and typed failures/cancellation retain both histories.
The focused fixture checks independent reciprocal current/motion and heat
balances, material-only substitutions, continued states and refusal.

The core fit stays at its separately held preparation temperature and biased
reference. This connection changes motion through the actual winding; it adds
no unsupported core thermal, hysteretic or acoustic law.

### Complete coolant on planar induction fields

All four original planar induction variants consume this same complete
`HeatAssembly` through `fs_couple::circuit_enthalpy::eddy::thermal`.
Fixed linear and nonlinear magnetic laws retain their native field equations;
explicit thermal laws retain field entropy, internal energy, parameter work
and reversible heat. Each proposed thermal endpoint resolves original
conductor sources and finite-fluid transport together. Passive coolant cells
retain their original caloric storage outside the magnetic mesh.

Native Ohmic heat enters its solid receivers once, nonlinear numerical loss
stays separate, and closed-enclosure roundoff is excluded from external heat.
The full fluid/heat payload and typed provider failures survive. Prepared
circuit bindings and resolved electrical work without a field owner refuse.
Original source, geometry, field and combined-capacity limits remain active;
rejected candidates preserve both histories.

### Complete coolant on three-dimensional induction fields

The original linear, fixed-law nonlinear and temperature-dependent nonlinear
`fs_couple::circuit_enthalpy::spatial_induction` drivers accept this same
complete `HeatAssembly` through `heat_assembly_candidate`,
`nonlinear_heat_assembly_candidate` and
`nonlinear_thermal_heat_assembly_candidate`. Finite coolant remains in its
native caloric cell. Each field/thermal trial resolves transport at the same
proposed temperature and retains the full accepted fluid sample with the
original magnetic field history.

Local Ohmic heat enters its original whole-volume receivers once. Affine-mu(T)
and explicitly bound separable nonlinear temperature laws retain their
original magnetic entropy, internal energy, parameter work and reversible
heat. Backward-Euler loss remains separate from physical heat. Native source,
charge, volume, field-domain, capacity and radiation-energy gates remain.
Graph circuit work cannot be discarded through this thermal provider, and
typed source refusals or cancellation publish no partial field or heat state.

The finite-links spatial-induction fixture uses an independent Whitney-mode
and finite-fluid balance for source-only viscosity, heat-capacity and
resistivity changes, continued intervals and transactional refusal. Its
nonlinear checks use fixed- and affine-mu(T) linear-material limits and
independently account for numerical energy export. No moving geometry,
arbitrary thermal magnetic law or empirical device validation is implied.

### Complete material heat on the native explicit timeline

`SourceConvectionHeatNetwork::observe_with_links` composes the complete prepared
heat assembly with an existing `DepositionObserver`. Original thermal cells
retain their mass, caloric laws, scheduled events and enthalpy history. One
starting-state evaluation supplies material conduction, contact/layer resistance,
closed radiation, prescribed-reservoir powers and finite-fluid link conductances.
Finite-fluid transfer stays reciprocal internal exchange on native links;
prescribed convection contributes external heat once. No additional thermal
advance or circuit is introduced.

The same native timestep owner combines its link restriction with the summed
per-cell prescribed-convection hA and frozen-emissivity radiation bounds.
The observer retains the original source actually applied at the start of the
last accepted substep separately from a fresh accepted-endpoint sample. Events
invalidate pending and endpoint data before the next complete active-set
evaluation. Splitting an advance at an actual event boundary preserves native
histories and the generated-coefficient sequence.

Thermal-only admission is shared with `PreparedHeatNetwork::observe`: graph
resistors, electrical contacts and circuit-basis bindings refuse rather than
discard electrical work. Competing link owners, source-domain failures,
cancellation and work exhaustion discard the thermal/observer candidate.
Callbacks must be repeatable; cloning must deeply copy mutable physical history.
External effects or shared mutable callback state cannot be rolled back.
Signed prescribed heat and the bounded closed-radiation numerical net remain
separate observations of the original native energy account.

Integrated conduction layers can coexist with disjoint finite-fluid links via
`with_integrated_conduction_and_link_conductances`. Their original
identity/geometry checks remain; an isotropic integrated layer cannot replace
an explicitly directional conductivity binding.

The material-substitution regression in
`tests/material_cell_convection/film_source/heat_network/finite_links/observer.rs`
compares native trajectories against independently written flat-plate and
explicit heat-balance equations. It covers finite coolant storage, a prescribed
reservoir, a scheduled heat source, event continuation, source provenance,
ownership conflicts and candidate rollback. Nonzero enclosure radiation has
its own finite-cell regression in fs-couple. These remain uniform-cell,
coefficient-frozen explicit steps: refine the timestep for temporal accuracy.
Flow speed, pressure and composition keep their original authored/provider
ownership; this adapter adds no fluid advection, momentum or pressure-work law.

### Forced spheres on the original finite-cell heat owners

`SourceCellConvectionLaw::ForcedSphere` admits Whitaker convection in both
prescribed-reservoir cell boundaries and reciprocal finite-fluid links.
The primary property convention is explicitly `Reservoir`: the actual
fluid bulk supplies density, Cp, conductivity and viscosity for Re, Pr and
h. The same source resolver is additionally queried at the current wall
temperature; only that second fluid state's viscosity supplies the original
bulk/surface correction. Mean-film substitution is refused.

Both original material owners and their checked thermal views remain in the
result through the optional `surface` sample. A shared query admission checks
requested temperature and agreement between each owner's thermal and momentum
transport views. The caller still binds one immutable source law, pressure and
composition across these queries; a generic transport trait cannot prove that
identity. Source-domain, viscosity-ratio and Reynolds restrictions propagate.

Characteristic length is the physical sphere diameter, and exposed or linked
area cannot exceed its complete surface `pi D^2`. Contact coverage and
equal-and-opposite link heat retain the native thermal owners. Each represented
sphere/patch must satisfy the isolated, isothermal correlation assumptions;
this is not a local law for arbitrary nonisothermal curved surfaces. Neither
a stagnant-gas conductance nor a second thermal store is added.

The existing implicit heat/circuit residual and explicit source/link observers
therefore re-query both property temperatures at their own actual trials and
substeps. Four focused G1/G3/G4 regression sources cover independently
calculated liquid-sphere heat, original source identity and domain refusals,
implicit reciprocal energy and reversed endpoints, and native observed link
advancement. Rust syntax/format and scoped whitespace checks passed; these
new integration regressions were authored but not executed in this change.

### Variable-viscosity turbulent circular tubes

`evaluate_tube_from_material` retains separate original bulk and wall-temperature
fluid sources for the shared `SiederTateTurbulent` card. Bulk density, heat
capacity, conductivity and viscosity supply Re, Pr and the dimensional `Nu k/D`;
only wall viscosity enters its correction. The source-backed `Forced` cell,
finite-link and reservoir paths automatically query the SAME fluid at both
actual temperatures when this card is selected. `Reservoir` property sampling
is required; a film tuple, missing heated length or rectangular aspect ratio
refuses. Every original source/domain and caloric check remains active.

Native liquid/gas graph heat bindings use the same two snapshots at actual
step-average flow and recheck the initial/final rates and any zero crossing.
The authored circular diameter must agree with the original pipe area. Wall
density is not a second stored fluid state; the original bulk density still
owns inventory admission. Finite walls, source-integrated wall surfaces and
prescribed reservoirs keep their original heat owner, conservation and clock.
No diagnostic watts are added again as a source.

`finite_gas::tube::ForcedTubeGasSurface` carries that same source-pinned wall
viscosity into the original finite-body piston and fixed-volume wall-feedback
solves. Its circular diameter and full heated length determine L/D and bound the
inside-wall patch area. `GasHeatLaw::ForcedTubeConvection` uses the same resolver
for reacting and transported-composition networks; the actual trial composition
and selected source claims remain authoritative. Both heat directions and body
latent heat retain the original reciprocal energy owner. Explicit homogeneous
Peng--Robinson branches retain real bulk density and Cp. This is the existing
engineering closure with bulk density, Cp and k frozen for the correlation,
not empirical gas validation or support for strong property variation,
compressibility, critical-region or fluid phase-change heat transfer.

The correlation is restricted to its declared smooth, straight, fully developed,
single-phase Newtonian regime. Geometry, pressure/composition identity and the
frozen-property approximation remain caller obligations; generic borrowed
snapshots cannot prove common material identity. Original source re-evaluation
and native energy closure do not establish empirical accuracy, turbulence
resolution or temporal order. Focused regressions cover both viscosity-ratio
directions, original source queries, reciprocal enthalpy feedback and refusal.

### Material-backed bare tube-bank crossflow

`medium::convection::evaluate_tube_bank_from_material` connects the original
`fs-convection::tube_bank` mean to gas/liquid transport and conduction Robin
boundaries. The caller supplies the same fluid at mean-bulk and wall
temperatures plus explicit inline/staggered geometry and approach speed.
Each source is sampled once and retained. Bulk density/viscosity give Reynolds,
bulk Cp/viscosity/conductivity give Pr and dimensional h, and the wall-fluid
snapshot supplies only Pr_surface. The shared geometry computes the actual
maximum-gap speed, including diagonal constrictions; source/domain refusals
and unquantified discrepancy remain unchanged. Robin lowering rereads neither
material and uses the retained bulk conductivity and outside diameter.

The generic transport seam cannot prove common composition/pressure identity.
The caller retains that obligation and the bank mean's effective bulk reference;
this snapshot API supplies no finite-stream heat balance, pressure loss, short-pin or
finned-surface correction. Its focused material regression checks distinct
bulk/wall snapshots against independent formula arithmetic, source identity,
gap Reynolds, stagnation refusal and the actual Robin coefficient.

`tube_bank::stream::LiquidTubeBankStream` adds a finite single-phase liquid
stream around one uniformly isothermal bank. It binds the original
`LiquidTemperatureLaw` and `IsobaricHeatCapacityLaw`, requiring the same card,
selected Cp claim, typed temperature axis and complete fixed source context.
Source scalar/temperature curves must be continuous and supported over the
inlet-to-wall interval; no material extrapolation or new caloric chart occurs.

The original source-caloric `linear_heat_exchange` owns the complete Cp(T)
integral and inverse. A bounded scalar bracket refreshes actual liquid transport
at the arithmetic inlet/outlet mean and surface Pr at the held wall. One mass
flux supplies the representative approach speed at mean density, so Re uses
mass flux and the actual mean viscosity without an extra inlet/mean density
factor. The actual inlet velocity is retained separately. Tube count, original
transverse pitch and axial span determine frontal area and lateral heated area;
half-pitch side margins or transverse repetition exclude inferred bypass area.

One retained relative-enthalpy increment supplies equal and opposite wall/fluid
watts. Original Cp at the represented outlet and the effective caloric Robin
mean provide two independent watt residuals; explicit absolute/relative gates
apply before any result returns. The result retains the original caloric owner,
complete inlet/mean/wall source states, correlation evidence, iteration residual
and the same actual conduction Robin handoff. Refusal or cancellation returns
no partial state and changes no source owner.

This is a steady, mean-property, isothermal-wall closure. It does not establish
empirical accuracy under strong property variation, tube-end insignificance,
uniqueness of nonlinear mean-property solutions, pressure loss, wall/fluid
storage, fins, phase change or a second stream. Zero flow and out-of-domain
forced Re refuse instead of manufacturing a stagnant-flow coefficient. The
existing `CaloricHeatExchanger` remains the separate two-stream owner. Focused
regressions cover independent constant-Cp heating/cooling/equilibrium, a
manufactured variable-Cp/transport endpoint, conserved mass-flux Reynolds,
parallel-area scaling, caloric datum changes and source/work/cancellation
refusals. Native execution status is reported separately from authored tests.

### Reciprocal material convection on finite thermal links

`source::links::evaluate_source_link_convection` resolves the original fluid
material at the actual fluid-cell or wall/fluid film temperature and supplies
`G = h A` to an existing thermal link. Both endpoint temperatures come from
their original caloric laws. The native link applies equal-and-opposite heat,
retains its interface coverage, and participates in the existing thermal
conductance/timestep owner. No external source power, duplicate fluid storage,
mass update or separate time integrator is introduced.

One link has one conductance owner; repeated bindings, foreign fluid endpoints,
nonzero authored conductance, competing overrides, unsupported material states
and cancellation refuse. Distinct links may share a finite fluid cell. The
source resolver must retain the same fluid material/conditions as its caloric
law; a generic transport trait cannot establish that identity. Fixed geometry,
uniform-cell validity, prescribed relative flow and correlation/source support
remain explicit. Focused regressions derive liquid storage and transport from
one material card, independently solve a closed wall/liquid backward-Euler
balance, change viscosity data, reverse link orientation, and check zero external
heat plus transactional source/ownership refusals.

`source::links::observer::SourceLinkConvectionObserver` installs these
conductances on the original explicit thermal timeline. Both the native pair
exchange and its combined conductive CFL consume the same values. Reports retain
the actual starting-state material evaluations used by the accepted substep and
fresh endpoint evaluations separately. Source-support failure at that endpoint
or cancellation discards both finite-cell histories and the candidate observer.
Competing conductance owners refuse through the existing observer composition.

The adapter preserves explicit integration and coefficient-freezing refinement
requirements. Focused tests change only liquid Cp source data, rebuilding its
original caloric law from the same card, and compare both finite temperatures
against the analytic constant-property solution across three timestep sizes.
They also check native CFL budget enforcement, zero external heat, distinct
applied/endpoint sources, late support failure and cancellation.

### Enclosed horizontal, inclined and vertical air layers between finite walls

`source::enclosure::evaluate_source_enclosed_gas_layers` supplies total gas-gap
conductance to an original reciprocal thermal link. Both endpoints are finite
isothermal walls; the gas has negligible storage in this quasi-steady model.
The declared lower wall determines signed buoyancy independently of endpoint
ordering; at 90 degrees it names the reference wall for signed heat. Actual
caloric wall temperatures determine the full contrast and the
mean temperature used to query one original gas material, retained with its
thermal snapshot and original correlation evidence on every trial.

At `inclination_degrees=0`, the shared `HollandsHorizontalAirLayer` card owns the onset and gas conduction
limit. Stable heating from above, subcritical heating from below, zero gravity
and exact equilibrium retain `Nu=1`; an unstable admitted layer adds buoyant
transfer. This is TOTAL gas heat: no second Fourier gas path or generated pair
of source rates is added. The native link alone debits and credits the same
heat and includes its conductance in the original timestep bound.

Positive wall-plane angles through 60 degrees use the shared
`HollandsInclinedAirLayer` card with the FULL gravity-based Rayleigh number.
Its normal driving and along-slope sine correction remain inside the original
correlation owner. Merely projecting gravity into the horizontal relation is
not the inclined model. A negative signed Rayleigh number at these inclinations
refuses; the horizontal stable conduction limit cannot establish an inclined
heated-from-above model. The actual declared angle is retained with each
source sample and correlation evaluation. The original observer and complete
heat assembly consume this same evaluator, including during nonlinear trials.

Exactly 90 degrees selects `WrightVerticalAirLayer` using the magnitude of the
full gap-based Rayleigh number. Either wall may be hotter; the original link
retains the signed temperature difference and reciprocal heat. `plate_spans`
then means `[whole width, whole vertical height]`, with height/gap in `[5,110]`
and width/gap above 10 as a conservative lateral-extent gate. The shared card
admits `Ra` in `[0,1e6]` and `Pr` in `[0.70,0.72]`; equilibrium retains `Nu=1`.
The second span is the physical vertical height, never the smaller span or a
mesh-cell length. Angles strictly between 60 and 90 degrees remain unsupported.

Bindings declare the normal gap, both complete rectangular plate spans,
full gravity, wall-plane inclination in degrees and an explicit maximum
`|beta delta_T|` in `[0,1)`. The native link
must cover that full rectangle on both walls within floating-point roundoff and have zero
authored conductance. Orientation-specific geometry, air-only source support and all
Rayleigh/Prandtl limits remain with the original card. A matching Pr is not
evidence that a non-air fluid is supported. Wall temperatures must represent
uniform isothermal faces; this adapter does not replace resolved wall gradients,
gas thermal storage, unsupported inclinations, tilted heated-from-above
circulation, sidewall bridges or radiation. The caller must bind the actual
wall orientation; the link's endpoint ordering is not an inclination model.

`source::enclosure::observer::SourceEnclosedGasLayerObserver` reevaluates the
gas at every actual explicit substep and accepted endpoint, retaining applied
and endpoint samples separately. The same evaluator can supply a native
backward-Euler or coupled heat trial's conductance map. Competing link ownership,
geometry/source-domain refusal and cancellation leave no partial map or state;
the original observer clone contract governs rollback of callback history.
Pressure and composition remain with the original gas resolver, and numerical
and physical refinement are still required for the declared approximation.

`source::enclosure::heat::SourceEnclosedGasHeatAssembly` carries those same gas
layers through the complete native `HeatAssembly` interface. Evaluate the gap
bindings at the accepted state, pass their conductance map to
`MaterialHeatNetwork::prepare_with_link_conductances`, then retain that prepared
owner and the same resolver in the assembly. Every trial re-evaluates the gas
from both actual wall temperatures before the prepared network resolves its
original material resistance, contacts, conduction and closed radiation. The
accepted `EnclosedGasHeatSample` retains the complete `HeatNetworkSample` and
all gas-source/correlation samples. No gas store, extra source pair or second
integration is introduced.

The original link partition refuses missing or overlapping Fourier ownership,
and a gas-gap link cannot also carry material contact resistance. Current
contact pressures and disjoint path metrics reach their original owners; a
metric affecting either gap wall refuses without explicitly updated gap,
plate spans and orientation. Native link conductances already enter the
thermal timestep bound. Additional explicit-source bounds delegate to the
original radiation owner, so the gas contribution is not counted twice.

Wrap the assembly in `SourceConvectionHeatLayer` to add outside prescribed
reservoirs on the remaining exposed wall area. The original circuit/field
residuals and `observe_heat_assembly` consume the complete assembly and retain
its nested samples. Gas-gap transfer and closed radiation remain internal to
the two-wall system; native circuit loss heats its declared material cells,
and only prescribed reservoirs or other authored supplies contribute external
heat. Existing source, work-budget, energy and cancellation gates retain their
original ownership and transaction boundaries.

Vertical-layer tests reconstruct Wright's expression independently, exercise
both reference walls and heat directions, distinguish vertical height from
width, and retain film-property refresh. The existing native backward-Euler
and coupled circuit/radiation tests also include vertical gaps, checking both
wall storage changes and the complete energy balance. Manufactured source packs
exercise the material path; they do not establish empirical fluid validation.

### Shared thermochemical medium in hydraulic and heat paths

`FrozenFlowMedium::bind_fan_bank` connects the same original reference/inlet
gas owners to the existing `FanBank`. Static pressure and its slope gain the
explicit density-affinity factor; volume-flow coordinates, source samples,
speed law and stall limits remain unchanged. The shared pressure factor is
used by the certified enclosure root, arbitrary graph operating point,
series/parallel composition, and their local speed/curve derivatives. Retained
material snapshots participate in operating/composite identities, and the
adapter keeps the original material/source-pack owners. A detached `FanBank`
retains numerical state and applicability, not new material-source authority.

The caller declares positive density and rotor-Reynolds ratio bounds and a
static-pressure fraction in `(0, 0.1]`. Admission checks coherent gas states and
the full tolerance-inclusive pressure envelope; liquids, repeated binding and
unsupported similarity refuse. Local derivatives also require interior
Reynolds and pressure support. Material fan composition requires the same
admitted inlet state on every member, preserves it through nested composition,
and checks the combined pressure against the strictest common limit. It cannot
mix unbound inputs or apply density twice. The evolving compartment fan reuses density
admission and rejects a frozen binding to avoid scaling twice. This is a
frozen, incompressible static-pressure similarity model: it supplies neither
total-pressure conversion, fan efficiency, material uncertainty propagation,
Mach/stall/chemical-compatibility validation nor measured exotic-gas accuracy.
`tests/material_fan.rs` checks source-backed pressure/composition/temperature
substitution through fan, hydraulic and heat calculations, analytic graph
roots and derivatives, composition, and refusals.

`graph::thermal::medium::FrozenFlowMedium::ThermochemicalGas` retains the
original `ResolvedThermochemicalGasState`, including species, caloric, EOS and
transport sources. Its single transport snapshot supplies hydraulic density,
Re/Pr, convection, Robin boundaries and the series path's conserved heat-capacity
rate. Source pressure or composition changes therefore reach the actual operating
point and heat march without rebuilding a legacy gas specification. These paths
remain frozen-property, constant-density sensible-heat approximations; they do
not advance chemistry, gas energy or phase, and they retain correlation-domain
refusals. `tests/material_medium.rs` checks solved pressure-driven flow, serial
heat balance, data-only pressure/composition changes and source-object retention.

### Material-driven hydrostatic pressure and natural circulation

`series::hydraulic::pressure::SeriesGravity` compiles a connected path and an
arbitrary uniform gravity vector. `SeriesPressureFlow::solve_with_gravity`
balances `sum(dp_loss) = dp_imposed + sum(rho_i g.dot(dx_i))` with the same
single source snapshot already used by each hydraulic region. Actual gas or
liquid density controls the head; no ideal-gas expansion, air constant, or
reference volume flow is substituted. Signed flow, zero imposed pressure,
opposing pressure, and zero gravity use the original resistance inverse.
Density-reference summation with the path's endpoint potential makes a
constant-density closed loop exactly quiescent and preserves small density
contrasts. The original convection handoff consumes the resulting local flow.

`LiquidHeatPath::solve_pressure_gravity_driven` recomputes bulk source density
and viscosity on every original thermal-path trial before balancing pressure.
The same bounded root thus selects the flow, downstream enthalpy and wall heat
with no pump pressure when gravity supplies a positive drive. Film temperatures
remain the convection property convention; they never replace bulk density in
the hydrostatic integral. This overload admits authored hydraulic losses, not
the separate pressure-dependent porous-face model. The positive mass bracket
selects one orientation; it is not a uniqueness or stability proof. The
directional circulation entry point below constructs the reverse traversal
without replacing the original material source or closing inlet reservoir.

Results retain local heads and signed conservative gravitational work
`mdot*g.dot(last-first)` separately from hydraulic dissipation. Net gravity
work on a closed loop is zero even when density contrast drives circulation;
neither head nor work is silently deposited as heat. Frozen-state paths do not
evolve an EOS. The prescribed-inlet solve remains a regional bulk-density, single-phase,
quasi-steady approximation with a prescribed inlet reservoir and its existing
explicit hydraulic-heating omission budget. It does not close a finite loop's
thermal inventory, compressible acceleration, pressure-dependent phase changes,
or global mechanical/thermal energy. G1/G3 regressions cover analytic flow and
enthalpy closure, source temperature/pressure changes, reversal, rigid rotation,
translation and exact uniform-loop equilibrium; G4 checks bounded refusal.

`LiquidHeatPath::reservoir` and `advance_reservoir_with_gravity` close a finite
well-mixed coolant inventory around that same gravity operating point. The
existing backward-Euler reservoir owner solves
`M*(h_new-h_old) = dt*mdot*(h_return-h_supply) + external_heat` on the original
caloric chart. Each inlet trial recomputes source transport, regional bulk
density, gravity head, hydraulic losses and the actual positive circulation;
the final circuit retains exactly the accepted inventory state as its inlet.
The template inlet remains immutable. Original source identity, mass and time
admission are shared with the existing pump-reservoir path.

The path must return to its exact initial vertex, not merely the same elevation.
This makes net conservative gravitational work zero and excludes an unmodelled
reservoir elevation transfer. `GravityReservoirControls` shares actual flow
evaluations across every inlet/pressure trial and enforces both the original
hydraulic watt allowance and the reservoir step-joule allowance. Source failures,
missing positive roots, exhausted budgets and cancellation leave both inventory
state and time unchanged. Only the reservoir stores heat; prescribed walls are
external thermal sources, the loop remains quasi-steady, and no pipe residence
time, inertia, continuous passage through zero flow or zero-flow conduction is
supplied. Tests compare repeated inventory steps against independent analytic heat/pressure
balances and exercise work accounting, closed geometry and atomic refusal.

`LiquidHeatPath::solve_gravity_circulation` supports explicit forward/reverse
selection and `LiquidFlowSelection::OnlyBracketed`. Reverse traversal rebuilds
the original region, wall, resistance and vertex order together and negates
the imposed pressure; it recomputes the complete source-resolved thermal path.
`DirectedLiquidSolution` retains receipts in actual traversal order and exposes
mass flow, pressure drop, head and residual with signs in the authored direction.
Heat and dissipation retain their physical signs. `OnlyBracketed` evaluates both
declared orientation brackets under one actual flow-evaluation allowance and
accepts exactly one bracketed candidate. Two candidates are an explicit
ambiguity, no candidates refuse, and source/correlation failures never become
evidence of an absent root. This is not global root uniqueness or stability.
`advance_reservoir_circulation` uses the SAME inventory integrator and atomic
commit for either direction, with work counted across both orientations and all
inlet trials. Accepted steps can switch orientation as imposed conditions change;
the crossing through zero and pipe momentum/storage remain unmodelled. G1/G3/G4
tests cover independent directional heat/pressure balances, branch ambiguity,
unsupported-source refusal, shared work and energy-preserving inventory reversal.

`FrozenFlowMedium::ThermochemicalRealGas` retains the original nonideal source
and uses its PR density and complete Cp through the same hydraulic losses,
channel heat-capacity rate and convection owners. Every frozen-medium variant
also delegates the common `ThermalFluidTransport` view.
`evaluate_vertical_plate_from_thermal_material` obtains reservoir temperature
and isobaric expansion from that same source, checks thermal/momentum tuple
agreement, and executes the existing Churchill-Chu correlation and Robin seam.
Unknown expansion and zero/out-of-range buoyancy refuse instead of acquiring
an air value. Coefficients stay frozen at the reservoir state; no film-state
reconstruction, phase evolution or compressible channel law is supplied.
The separate fan-similarity owner still requires ideal `GasState` snapshots:
nonideal inlet and reference bindings both refuse instead of substituting the
retained ideal-reference source for the actual real-gas medium.
Focused tests cover source-pack pressure/composition changes through native
channel heat marching and natural convection, the ideal-gas explicit-path
limit, inconsistent snapshots and unknown-versus-zero expansion.

`medium::steady_gas::solve_material_gas_transport` closes the native
`SteadyAxialGas` stationary-species/NASA9 enthalpy solve against convection
evaluated from each retained local source state. The same downstream mixed-cell
state supplies velocity, Re/Pr, heat transfer and accepted hydraulic loss;
composition injections retain their native mass and enthalpy transport. The
returned `AxialGasField` remains usable by existing gas-field consumers.

Acceptance bounds the native power residual plus the actual wall-law defect,
with bounded feedback iterations and a total native-substep cap. Hydraulic
Reynolds and pressure-drop admission use the converged thermal state, so a
nonphysical zero-wall initial guess cannot reject an admitted final flow.
Cancellation and any refusal publish no partial field and do not mutate the
input. `tests/steady_material_gas.rs` checks independent heat/pressure balances,
temperature-dependent transport, species injection, pressure changes,
final-state hydraulic admission and transactional failure/cancellation.

`solve_material_gas_network` uses that same material-convection fixed point for
the native branched and inlet-fed recirculating `SteadyGasNetwork`. There is one
feedback implementation for the axial and graph APIs; each keeps its original
thermal result type. Each cell's full accepted source state supplies its heat
transfer, hydraulic response and acoustic medium. Composition and temperature
changes therefore update the wall law inside the steady balance, rather than
leaving an independently authored coefficient frozen. Correlation evidence and
material receipts remain bound to the returned local states.

Cells must represent the declared bulk-flow sections and correlations; junction
or nozzle geometry is not inferred. Per-cell hydraulic losses are observations
of prescribed flows, not a graph pressure solution. Their sum is a conservative
common-pressure admission bound; on an axial series path it is the hydraulic
loss contribution, excluding density-driven acceleration. These prescribed-flow
APIs add no pump or flow-feasibility claim. The network regressions reconstruct cell and external enthalpy
balances, check recirculating species, exercise material/pressure substitutions,
and compare the original axial API with the same graph and wall laws.

`medium::convection::finite_gas::ForcedGasSurface::isochoric_body_candidate`
connects the same material/correlation evaluator to the generic finite-gas,
uniform-body thermal callback in `fs-couple`. Each coupled trial resolves h from
the actual candidate gas endpoint; accepted transport and correlation sources
remain beside the same gas state used for energy and acoustic observations.
The original body solver retains conductivity, radiation, Biot admission and
enthalpy, and the single exchanged heat cancels in the combined energy balance.
The same isochoric adapter admits both ideal and Peng--Robinson parcels through
the shared material storage interface. Real EOS density and Cp supply Reynolds,
Prandtl and the live coefficient; the accepted parcel retains its real pressure,
caloric storage, entropy and acoustic state. No ideal-reference density enters
that transfer. `tests/finite_gas_convection/real_gas.rs` covers dense-gas heating
and cooling, reciprocal energy/entropy accounting, wall feedback, formation
reference invariance and transactional source/cancellation refusals. Film
queries use the explicitly selected ideal or homogeneous Peng--Robinson source.
A positive patch must cover the whole body surface; zero area explicitly
insulates it. Bulk-gas property evaluation and prescribed relative flow retain
their correlation limits. Dittus-Boelter selects its heating/cooling branch
from the proposed signed heat inside the same residual. Its neutral zero-heat
seed uses the smaller branch coefficient, preventing an authored direction
from falsely rejecting the physical branch at the Biot gate. Constant-heat-flux
duct correlations refuse on this isothermal-body path.

`ForcedGasSurface::finite_body_candidate` also accepts the original real-gas
piston through its shared conservative thermal/motion candidate. The actual
Peng--Robinson midpoint supplies density and complete Cp for Reynolds, Prandtl
and the forced-convection coefficient, while its temperature selects heating
or cooling. That same midpoint supplies pressure work in the mechanical solve;
the finite-body heat debit and gas heat credit remain one reciprocal exchange.
The returned sample retains the concrete gas material, transport and original
correlation evidence. The moving-gas regressions exercise both heat directions,
independent total energy, retained source views and atomic refusal. Relative
speed remains prescribed; no resolved flow or fan-energy account is inferred.
The explicit constant-h API remains available. `tests/finite_gas_convection.rs`
checks endpoint coefficient re-evaluation, transport-only data substitution,
energy closure, fixed-h identity and transactional refusals. Additional tests in
`tests/finite_gas_convection/uniform.rs` exercise pressure/composition changes,
latent heat, both heat directions and late source-callback cancellation.
No spatial flow,
film-temperature solve or experimental fidelity claim is introduced.

`solve_material_gas_transport` is a first-order steady reduced model with prescribed positive mass flow,
supplied wall temperatures, downstream mixed-cell bulk properties and a common
static EOS pressure. Total pressure loss must fit the explicit small-drop
limit. It does not solve fan mass flow or axial pressure, advance reactions or
phase, resolve boundary layers, or claim experimental accuracy. Hydraulic
`dp*Q` is reported separately and neglected in the thermal equation; it is not
silently added to wall heat or enthalpy.

### Source-preserving property temperatures in finite gas/body coupling

`finite_gas::film::FilmGasSurface` evaluates external forced or natural
convection at the mean of the actual bulk gas and wall temperatures.
`FilmGasHeatLaw::BulkForced` instead selects the actual bulk temperature for
the existing circular/rectangular isothermal duct cards, including Hausen,
Dittus-Boelter and Gnielinski. The original correlation's geometry, flow regime
and modest-property-variation assumptions remain necessary. Bulk sampling is
the convention stated after Eq. (6) of [Hata and Noda's primary study](https://www.jstage.jst.go.jp/article/jsmeicone/2007.15/0/2007.15__ICONE1510_15/_pdf/-char/en);
it does not supply a missing wall-viscosity correction.

The existing finite-gas/body/piston residual owns both heat debits, gas pressure
work and mechanical motion. The property state supplies coefficients only;
the original bulk gas retains mass, energy, composition and acoustic state.
Rayleigh forcing and signed exchange retain the full wall-to-bulk temperature
difference. Duct cards passed to mean-film `Forced` refuse with guidance to use
`BulkForced`; external plate/cylinder cards remain in `Forced`. Active
constant-heat-flux duct cards refuse in the shared `ForcedGasSurface` admission,
including piston, isochoric and network paths, because the body is isothermal.
Zero-area insulation remains explicit.

The bulk-property path still re-queries the source bank at every actual gas
trial. It cannot bypass composition, original source or real-gas branch checks
by cloning the supplied snapshot. Regression source checks variable properties
and actual heat direction against three original duct equations, wrong-source
refusal at equilibrium, external-film preservation, reciprocal finite-body
energy and incompatible heat-flux laws. These new Rust tests were not executed
in the implementation session because the execution environment was unavailable.

`ThermochemicalGasModel::evaluate_temperature_from_snapshot` reconstructs
the actual active species and declared composition basis from the snapshot,
pins its original transport claims, and admits the original caloric sources
before changing property temperature at unchanged bulk pressure. A source bank
may have different initial fractions or additional inactive species, but cannot
substitute material cards, caloric laws or source packs. Source-domain failures
remain failures. This path retains the original dilute ideal-gas, correlation,
uniform-cell and temporal-refinement limits; it adds no finite film store,
boundary-layer solver or chemical evolution. Focused finite-gas regressions
exercise mixed and pure species, mass-basis preservation, independent film
Rayleigh/heat algebra, actual piston/body energy closure and atomic refusal.

`PengRobinsonFilmSource` extends that same surface and reciprocal heat/motion
owner to real gas with an explicit Vapor or Liquid branch. The material owner
pins the original PR critical properties, both interaction-matrix axes and
actual active source composition; ambiguous Stable selection or a bulk-root
mismatch refuses. External forced, oriented natural and aiding mixed
correlations use the actual film density and Cp; `BulkForced` instead re-queries
their bulk counterparts under the same source and branch checks. Buoyancy uses its source-owned isobaric expansion
derivative, including real-gas departure from `1/T`. The full wall/bulk contrast
remains in Rayleigh forcing and signed heat. A shared immutable PR law avoids
rebuilding its vectors and interaction matrix during ordinary thermal trials.
Focused regressions independently reconstruct the dimensionless groups and
correlation results, close gas/body/mechanical energy in both heat directions,
and preserve live storage on phase/source/cancellation refusal.

`examples/real_gas_heated_chamber.rs` joins the original real-gas nozzle,
finite reservoir, finite thermal body, PR film convection and moving piston.
Its CSV trajectory reports changing densities and inventories, heat, work,
mass defect and a common-reference whole-system energy balance. Synthetic
sources, homogeneous branches, an isothermal solid body, prescribed relative
flow and first-order splitting remain explicit; timestep refinement is needed.

### Finite-gas streamwise plate sections

`FilmGasHeatLaw::ForcedPlateSection(ForcedPlateGasSurface)` applies the shared
laminar or leading-edge-corrected turbulent strip relation to the original
finite-body/piston, fixed-volume wall-feedback and coupled heat-network paths.
`start` and `end` are physical downstream distances from one leading edge;
a numerical thermal patch does not restart the boundary layer. The caller
retains the real patch area and must rebind distances after flow reversal.
Both positive endpoints must remain inside the selected original card's
domain; unsupported geometry, cross-regime strips and underflowed positive
starts refuse.

Each trial queries the existing source owner once at the actual wall/bulk mean
and bulk pressure, preserving actual composition and material claims. Bulk gas
still owns energy, mass and acoustic state; the film is a property evaluation.
The original thermal owner applies one signed heat exchange and retains latent
heat, pressure work, area/Biot gates, cancellation and atomic publication.
Zero area remains explicit insulation without a film query.

`GasConvectionEvaluation::effective_nusselt()` returns the strip value. Its
`nusselt` field remains the original downstream mean; `section` retains the
upstream mean, physical distances and actual strip Nu. Moving these evaluations
into the coupled report avoids recomputation or metadata clones. Positive-start
coefficient evidence keeps an explicitly unbounded model discrepancy rather
than inheriting the mean card's uncertainty as a local error estimate.
At common isothermal wall/film states, the area-weighted strip heat equals
whole-plate heat. Local property refresh is a held-property approximation;
spatially non-isothermal boundary-layer history remains outside this model.

Focused `finite_gas_convection/plate_sections` regressions cover an independent
equal-half-plate integral, source/pressure/composition retention, both heat
directions, latent heat and gas/body/motion conservation, fixed-volume wall
feedback, insulation and transactional geometry/domain/cancellation refusal.

### Bulk and surface properties for forced spheres

`medium::convection::evaluate_sphere_from_material` retains both supplied fluid
owners and uses the shared Whitaker sphere card: bulk rho, mu, Cp and k supply
Re, Pr and h, while only viscosity comes from the fluid at surface temperature.
The helper also lowers the resulting coefficient to the original conduction
Robin row. The generic transport seam cannot prove that two supplied sources
describe the same fluid; its caller must preserve pressure and composition.

`finite_gas::sphere::ForcedSphereGasSurface` performs that second temperature
query through the existing source-pinned ideal/explicit-branch real-gas owner.
Its repeatable callback and finite-body adapter use the original reciprocal
heat/piston solver, with no additional gas parcel or thermal integrator. Gas
energy, body enthalpy (including latent heat), pressure work and motion retain
one signed heat exchange. Patch area cannot exceed either the sphere's pi D²
surface or its declared thermal body's exposed area. Zero area is insulation.
The complete thermal bracket must remain in the selected correlation domain.

Whitaker's source bounds, including mu_bulk/mu_surface in [1, 3.2], remain
strictly enforced; no film-property substitution, unit-ratio clamp, stagnant
extension, buoyancy blend, evaporation or extra conductive link is inferred.
Prescribed free-stream speed supplies no fan work or drag. Focused G0/G1/G3/G4
tests cover an independent liquid formula and Robin lowering, actual gas
source/pressure/composition retention, trial-dependent viscosity, melting and
gas/body/motion energy closure, deterministic replay, insulation and atomic
source/geometry/cancellation refusal. They do not establish experimental
accuracy or support outside the original isolated-sphere model.

### Wall-dependent fixed-volume gas/body coupling

`FilmGasSurface::isochoric_body_candidate` and
`ForcedSphereGasSurface::isochoric_body_candidate` connect these same laws to
the original fixed-volume `advance_with_wall_convection` transaction. Forced,
oriented natural and aiding mixed films use the current wall and bulk
temperatures; forced spheres retain bulk transport and separately query
wall-temperature viscosity. Property samples retain typed wall temperatures;
the original enthalpy owner retains the actual phase and latent fraction.
The existing wall solve bounds both temperature mismatch and boundary-work
mismatch before the original reciprocal heat, energy and entropy accounts
accept. Positive area must equal the body's complete exposed area, preserving
local Biot admission; zero area explicitly insulates. Constant-heat-flux duct
cards are refused for this isothermal-body boundary. The full heat and wall
brackets retain all original source, correlation and cancellation refusals.
Focused `isochoric_wall` regressions cover the coupled property states, full
Rayleigh contrast, energy closure, melting, partial-area/Biot refusal and
cancellation.

`FilmGasSurface::resolve_isochoric_coefficient` and its prescribed-mass-flow
counterpart expose that same source/coverage admission to
`FiniteGasBodiesEnvironment::advance_with_wall_convection`. Distinct bodies
can select different natural, forced or mixed surface laws while sharing one
finite gas endpoint. Each body retains its own queried and accepted wall/film
temperatures; all use the actual bulk pressure, composition and source claims.
Prescribed mass flow derives velocity from bulk density even when transport
properties use film density. The original thermal owners apply area and debit
reciprocal heat once. Focused regressions compare the single-body limit,
simultaneous hot/cold bodies, canonical input ordering, conservation, source
identity, mass-flux conversion, and atomic source/coverage refusal. This adds
no resolved flow, open mass/enthalpy flux or extra finite heat store.

### Natural convection in finite gas and coupled surface networks

`medium::convection::finite_gas::NaturalGasSurface` evaluates the shared
natural-convection geometry from the coupling owner's actual gas and body query
states. Its raw bulk-property entry retains vertical-plate behavior; active
horizontal plate/cylinder and sphere exchange require `FilmGasSurface`, which
evaluates the original ideal or explicit-branch PR gas at the mean film temperature.
The bulk vertical-plate path accepts ideal and Peng--Robinson gas states and
uses the existing `evaluate_vertical_plate_from_thermal_material` owner. Its
coherent thermal/transport view supplies actual temperature, density, Cp,
viscosity, conductivity and local isobaric volumetric expansion. The real EOS
derivative replaces the ideal-only beta = 1/T shortcut. Pressure and composition
therefore propagate through Ra, Pr, h and the original signed gas/body heat
residual. The same generic piston, phase-enthalpy and gas owners advance their
own energy, including latent heat and piston work.

`tests/finite_gas_convection/real_gas/natural.rs` reconstructs Rayleigh, Prandtl
and Churchill--Chu heat transfer from the actual source state, compares the
real derivative with the ideal shortcut, and exercises small hot/cold contrasts
through moving-piston energy closure. Exact equilibrium/isolation, correlation
and geometry refusals and late cancellation preserve both inputs. This local
derivative does not certify the bulk-property approximation over an arbitrary
temperature range or replace the existing correlation evidence and domain.

`finite_gas::network::exchange::GasHeatLaw::NaturalConvection` carries that
same law into the material electrical/contact/radiation/thermal network and its
chamber-flow step. Natural, forced and explicit conductive laws may occupy
disjoint exposed patches, with aggregate area admission and one energy owner.
Accepted samples retain the exact queried gas and wall states. Candidate solves,
atomic flow-step publication and cancellation use the existing network path.

`GasHeatLaw::{FilmConvection, ForcedSphereConvection}` carries the source-pinned
film and sphere laws through these same reacting, electrical, moving-wall and
conservative-transport assemblies. Film laws include forced, natural and aiding
mixed convection. Explicit property banks are queried at the actual trial's
pressure and composition: their initialization fractions never replace reacted
or transported gas. Film properties use the mean wall/bulk temperature; forced
spheres retain bulk Re/Pr/k and separately resolve wall-temperature viscosity.
The complete temperature contrast and original source domains remain binding.
Concrete property and source records survive in `GasHeatSample`; every patch
shares the original exposure budget, signed contact heat and atomic publication.
The network retains its existing ideal thermochemical gas storage owner.

`tests/finite_gas_convection/natural_network/source_laws.rs` and
`examples/material_convective_network/source_law_tests.rs` exercise simultaneous
sensible/latent exchange, changing transported and reacted composition, finite
circuit and wall-motion energy, source/domain refusal and rollback. They retain
admitted sphere Prandtl and viscosity-ratio ranges without clipping properties.

These are unconfined isothermal geometry correlations; they do not resolve
buoyant velocity or supply enclosure corrections. The caller must admit the
geometry, characteristic length and temperature contrast. Horizontal-plate
orientation and explicit `A/P` reach the same finite-gas heat residual, while
film properties leave bulk gas storage, pressure and composition ownership
unchanged. A focused horizontal-film regression closes the original gas/body
energy exchange.
The isolated sphere uses its authored diameter and the same film-source and
energy owners. Zero gravity or zero sourced expansion reaches its physical
steady conduction limit `Nu=2`, already included in `h*A`; no parallel
conductive path is added. The finite-gas zero-gravity exception applies only
to spheres. Focused sphere tests cover sourced Cp/k/density, both heat
directions, native finite-link storage and joint gas/body/piston energy.
Exact equilibrium and zero patch area have zero transfer without a fabricated
Ra=0 coefficient; every nonzero contrast retains the original correlation
domain. Tests in `tests/finite_gas_convection/natural.rs` and
`tests/finite_gas_convection/natural_network.rs` cover independent correlation
values, pressure/composition substitution, both heat directions, latent heat,
gas/body/piston energy closure, mixed surface laws, source identity, area
overlap and late-cancellation rollback. The flow-step regression uses zero
nozzle links and checks the original publication path.

### Mixed gas heat laws after combined species transport

`finite_gas::network::exchange::{transport_step_candidate, advance_with_transport}`
connects the same conductive, forced and natural patches to the shared
same-snapshot nozzle/Maxwell–Stefan transport owner. Both connection lists
consume one original species allowance and reconstruct each gas once. The
actual transported/reacted gas enters every heat/momentum trial, including its
composition, density, heat capacity, viscosity and thermal conductivity. No
ambient composition, coefficient cache or separate heat integrator is added.
One resolver supplies every mixed-law path; the original `fs-couple` owner
accounts for the complete interval and atomically publishes all histories.

Nozzle and molecular paths remain distinct physical declarations. Empty lists
disable either mechanism. Existing molar-average diffusion T/p restrictions,
withdrawal/entropy/source gates, convection domains and first-order splitting
remain binding. A forced patch still has prescribed relative speed; a nozzle
mass rate does not establish a resolved surface velocity or jet thrust.
`tests/finite_gas_convection/natural_network.rs` exercises nonzero nozzle and
diffusive transport with mixed convection, exact independently staged physical
endpoints, original source composition in accepted heat samples, signed cell
heat, full endpoint energy, diffusivity-only trajectory changes and complete
rollback after a late correlation refusal. Synthetic sources establish these
integration identities, not empirical gas or enclosure accuracy.

### Coupled total-enthalpy endpoint adjoints

`graph::thermal::coupled_transport::sensitivity::enthalpy` binds the concrete
`EnthalpyRobinResponse` to the existing mixed-air transport linearization. It
reuses the temperature response's exact port-name/order, HTC, area, effective
reference and branch heat-balance admission. The solid's complete physical
Jacobian includes reference-mass storage, nonlinear conductivity and any
ambient-radiation feedback; the air sees only original convective exchanges.
The endpoint pullback carries direct nodal specific-enthalpy seeds separately
from temperature/wall/heat/air objectives, preserving latent history where
`dT/dh=0`. It closes the implicit transpose interface equation using the
existing unrelaxed residual checks and stationary or fresh bounded IQN-ILS
updates, never differentiating primal solver iterations.

Returned previous-h, P1 source-density, assembled-load, supply-temperature and
convective `ln(HTC)` gradients include mixed-air feedback. Inner-column limits
apply per solid sweep and the interface budget bounds sweeps; cumulative
Krylov work is reported. An exhausted, invalid or cancelled solve publishes no
partial gradient and leaves the accepted endpoint reusable. Geometry, time
grid, reference masses, chart and material laws remain fixed; the original
pullbacks also hold flows fixed. This API does not add individual branch-flow,
fan operating-point, phase-chart, adaptive-grid or enclosure derivatives.

`pullback_flow_scale` and `pullback_flow_scale_iqn` additionally return the
derivative with respect to `ln(s)` when every signed branch and external
heat-capacity rate scales by the same positive `s`, at fixed proportions,
directions, topology, HTC, wetted areas and inlet temperatures. They first
close the complete enthalpy/interface transpose, including direct h carry and
ambient-radiation feedback, and reuse the existing air-only Euler contraction:
explicit degree-one air heat objectives minus the sum of **air** log-conductance
pullbacks with the accepted interface multiplier. The solid HTC derivatives
are not subtracted. One extra air reverse sweep adds no solid Krylov work and
preserves the original thermal gradient. Contraction failure/cancellation
returns no partial result. A consumer must independently justify that its
physical fan control follows this common-flow scaling model.
`tests/coupled_enthalpy.rs` checks real perturbed coupled solves, mixed
sensible/latent slopes, explicit latent carry, convection-only radiative
ports, feedback sensitivity and refusal/reuse behavior. Uniform-flow checks
perturb actual hydraulic pressure by `s^2` and resolve the coupled endpoint,
including ambient radiation, direct heat objectives and a latent state whose
temperature derivative is zero while its h-history flow derivative is nonzero.

### Coupled discrete thermal goals

`conjugate::goal::compare_discrete_goal` compares two supplied full solid
fields under one or more independent air paths sharing that solid. It binds
the same names, coefficient and wetted areas as the actual uniform Robin
traces, independently admits every branch through the existing temperature
and watt gates, and remarches the air at both fields. Its physical residual
is `R_c(T) = R_s(T) - B [r(T)-r_bound]`. The transpose uses the actual solid
material/contact Jacobian and an analytic reverse sweep of the exponential
air law. A fresh, bounded IQN history solves the reduced transpose equation;
neither the primal iteration history nor perturbed solver runs supply a
derivative. The complete coupled primal and dual residuals are recomputed
before publication, and inner Krylov work plus interface sweeps are retained.
The interface equation normalizes the goal by its largest nodal weight, so
its absolute stopping threshold is independent of the objective's units.

The output retains signed full-node residual contributions and the measured
linearization remainder. Consumers must retain that remainder, including
nonlinear material truncation and finite solve arithmetic. These are
observations between discrete fields, not continuum bounds. Flow, h, material
laws, mesh, fixed temperatures and contact resistance remain fixed; the API
adds no hydraulic, moving-interface, radiation or mixed-network derivative.
`tests/conjugate_goal.rs` uses actual P1/conjugate primals, an independent
effective-conductance adjoint oracle at high NTU, downstream serial heating,
multiple independent inlets, quadratic nonlinear-remainder scaling, objective
unit rescaling, and binding/reference/budget/cancellation refusals.

## Public types and semantics

- `FanCurve` owns typed `(VolumetricFlowRate, Pressure)` points, source identity,
  tolerance authority, a stall boundary, and the validity range for RPM affinity
  scaling.
- `FanBank` represents identical fans in series or parallel. At speed ratio
  `s`, flow scales as `s` and pressure as `s^2`; series pressure adds and
  parallel flow adds.
- `fan_law_source` exposes the scaling authority as structured
  `SourceProvenance`: AMCA Publication 201-02 (R2011), section 6.5.1, rather
  than the AMCA 210 / ASHRAE 51 laboratory test-method standard.
- `LossElement` uses a typed `LossResistance` in Pa/(m^3/s)^2. `LossNetwork`
  composes quadratic elements recursively in series and parallel. A loss
  element remains cardless unless its owner explicitly attaches a finite,
  nonempty `ValidityDomain`; `with_regime_validity` then binds the card name to
  `airflow.loss.<element>` and the card version to the retained source
  identifier. `LossNetwork::regime_audit_cards` and the enclosure wrapper
  expose those exact cards without silently deduplicating ambiguous names.
- `sharp_edged_orifice_loss` lowers a declared open area and typed air
  density to a `LossElement` through the sharp-edged thin-orifice model
  `R = ρ / (2 Cd² A²)`, with `Cd` the retained plateau constant
  `SHARP_EDGED_ORIFICE_CD` (0.61) cited through
  `ORIFICE_CD_SOURCE_CITATION`/`ORIFICE_CD_SOURCE_IDENTIFIER`. The declared
  `ORIFICE_RESISTANCE_UNCERTAINTY_REL` (0.15, engineering allowance) covers
  the sourced `Cd ∈ [0.60, 0.65]` large-opening spread, and every produced
  element carries a `loss_reynolds` validity card over
  `ORIFICE_PLATEAU_REYNOLDS` so the plateau assumption is audited after the
  solve. Non-finite or non-positive area or density refuses as
  `InvalidOrificeInput`.
- `EnclosureNetwork` cannot be constructed without a distinct
  `LeakageElement`; leakage is not an implicit constant.
- `solve_operating_point` returns an `OperatingPoint` with an interval-Newton
  unique-root bracket for the nominal model, weaker physical uncertainty
  estimates, per-terminal flows, and nominal leakage fraction. Flow bounds use
  the low-fan/high-resistance and high-fan/low-resistance corners. Pressure
  bounds separately solve the low-fan/low-resistance and
  high-fan/high-resistance corners, then evaluate each root with the same
  resistance that produced it.
- `OperatingPoint::correlation_handoff` converts one branch flow to typed mean
  velocity, computes Reynolds through role-tagged `fs-regime` dimensions, and
  produces `fs-convection::CorrelationInputs` without discarding evidence.
- The test-only plate-fin fixture composes that public seam through three
  Shah-London rectangular-channel cards and the real `fs-conduction` solve:
  one narrow simultaneously developing CWT table slice drives fin flanks,
  while the fully developed CWT/CHF limits provide companion checks and the
  CHF base-floor row.
  Its base plate and three fins are separate tetrahedral blocks joined by
  three named, card-backed bondlines; internal fin flanks and channel floors
  remain separately named Robin traces. Channel area, wetted perimeter,
  hydraulic diameter, aspect ratio, and length ratio all derive from the same
  mesh constants.
- `qoi::extract_thermal_qois` consumes a `ConductionSolution`, its exact
  `ConductionMesh`, an `OperatingPoint`, and one `ThermalQoiDeclarations`.
  It emits the five E05.10 families: deterministic maximum junction
  temperature, pressure drop, fan input power, surface mean/uniformity, and
  thermal margin. An absent requirement refuses; there is no default limit.
- `ThermalQoiDeclarations` groups everything the caller declares: the junction
  and surface regions, the cited fan-efficiency interval, the cited
  requirement, and the optional cited `DiscretizationReceipt`. Keeping the
  declarations in one borrow separates what the caller asserted from what the
  solvers produced, and keeps the entry point narrow as families are added.
- `ThermalRequirement` carries the *effective* limit together with a mandatory
  `SafetyFactorAuthority`. Following the project-layer doctrine, no universal
  multiply/divide rule is inferred: the declaring source owns how the finite
  factor (finite, at least one) was applied, and the supplied temperature is
  already post-factor. The factor and its policy bind into requirement
  identity and never enter the margin arithmetic, so this crate cannot derate
  a second time. Holding the effective limit fixed while changing the factor
  rebinds identity and leaves the margin value bit-identical; that equality is
  the executable meaning of single application.
- `DiscretizationReceipt` is the only route by which the junction maximum
  acquires a known `Discretization` term. It is a caller-declared
  refinement-ladder or goal-oriented half-width **in kelvin**; this crate
  neither computes it nor derives it from a solver residual. Admitting one
  populates the budget term only and leaves `NumericalKind::NoClaim` in place,
  because a declared ladder observation is not a certificate.
- The thermal margin is a functional of the junction maximum, not an
  independently constructed budget. Because `margin = limit - maximum` has unit
  sensitivity in the same kelvin, every junction-maximum term value transfers
  1:1 into the margin budget; provenance is rebound to the margin identity so
  the requirement and its factor still rebind it. The declared limit
  contributes no term of its own: it is exact by declaration, not measured.
  Propagation copies representations rather than composing them, so an
  `Unknown` source term can never be strengthened in transit.
- Material authority is visible in the evidence, not only in the digest. Every
  temperature QoI names its conductivity provenance as a model card
  (`fs-conduction:material-matdb-receipts` or `fs-conduction:material-declared`)
  and carries a provenance-specific `Parameters` gap reason. Both remain
  `Unknown` — neither authority supplies a propagation theorem — but "receipts
  are retained, propagation is not" and "the values are caller-declared with no
  receipt at all" are different claims and now read differently. Before this,
  `ProvenanceClass` reached only `solution_identity`, so two solves with
  different material authority hashed differently but were byte-identical to a
  reviewer.
- `extract_thermal_qois` validates the retained report as well as the field.
  `ConductionSolution` has public fields, and a resumed or deserialized
  solution is assembled rather than solved, so the seam refuses (a) any
  retained linear record with `converged_true == false`, naming the iteration,
  method and true relative residual — a QoI taken from an algebraically
  unconverged field is unsupported, not merely weaker — and (b) a report
  claiming `MatdbReceipts` while reporting zero receipts, which is
  self-contradictory. The production solver already fails closed on the first
  case (`ConductionError::LinearSolveFailed`), so that guard defends the seam
  rather than a live solver path.
- Every emitted `ThermalQoi<T>` carries both the existing `Evidence<T>` view and
  an `EngineeringUncertaintyBudget` with exactly one term for roundoff,
  solver/algebraic, discretization, geometry, parameters, boundary conditions,
  model form, and measurement. A missing propagation theorem/receipt is a
  named `Unknown`, never `Negligible` or zero.
- Every emitted `ThermalQoi<T>` also carries a `ThermalQoiKind` stating what its
  scalar MEANS beyond its dimension. Five of the seven emitted records are
  `Temperature` labelled `kelvin`, but the junction maximum and surface mean are
  `AbsoluteTemperature` while the spread, face-mean standard deviation, and
  thermal margin are `TemperatureDifference`. Unit and Rust type together cannot
  separate them, so a consumer must read this discriminant rather than infer
  semantics from a field name: only an interval may be negative, and only an
  absolute value has a meaningful Celsius reading. `fs_qty::semantic::QuantityKind`
  is the authority for this distinction; it is not reused directly because it
  carries no mechanical-power variant, and fan input power is mechanical.
- The kind is enforced, not merely recorded: a QoI whose kind forbids a negative
  scalar refuses if one is produced, so a below-absolute-zero junction maximum
  fails closed instead of travelling downstream wearing an ordinary unit label.
  A negative `thermal_margin` is admitted and is load-bearing — it is the design
  missing its requirement. The final-envelope audit demotes evidence but never
  rewrites a kind: what a scalar means is not a function of its authority.
- `ThermalQoiSet::audit_operating_envelope` is the mandatory final E05.10
  product gate. It requires exactly one consumed-card declaration for each of
  the seven emitted records, derives each incoming color from the actual
  `Evidence` receipt, runs one all-card/all-point `fs-regime` audit, applies
  each exact receipt to the matching eight-term budget, and returns the audit
  beside the updated values. Missing, duplicate, or foreign QoI declarations
  refuse; callers cannot supply the pre-audit color. Fully in-domain admission
  leaves the QoI set byte-for-byte unchanged, while any partial/out-of-domain
  envelope makes the affected model-form term explicitly Unknown under the
  exact receipt identity. Overrides remain acknowledgements only.
- `FanCurve::model_card` exposes the exact fan-law card already consumed by
  airflow evidence so product admission can audit the retained base-curve flow
  and speed-ratio axes without synthesizing a parallel card definition.
- The E05.10 battery also consumes the public Dittus-Boelter correlation card
  alongside that fan card. Its adversarial variant changes only Reynolds
  number, proving the correlation receipt is the sole demotion cause while the
  actual fan card remains in-domain.
- `ThermalQoiSet::audit_operating_envelope_with_cards` accepts owner-neutral
  `RegimeAuditCard` identity/validity projections for matdb and other native
  card schemas. The existing evidence-`ModelCard` method is an exact wrapper;
  parity is executable, so widening the card source cannot change receipts or
  budget demotions.
- The operating-point pressure/flow envelope populates the conditional
  boundary-condition term for pressure and power. The cited total-efficiency
  interval populates the fan-power parameter term. Both remain accompanied by
  an `Unknown` model-form term, so finite conditional bands cannot upgrade the
  synthetic/quadratic airflow model to validated product authority.
- Surface mean uses exact P1 triangle integration (`area * vertex mean`).
  Spread is the selected surface-vertex range. The reported standard deviation
  is the area-weighted dispersion of face means; it is explicitly not an exact
  integral of the pointwise squared P1 field.

## Invariants

1. Fan flow points increase strictly and pressure never increases.
2. Fan-law scaling refuses outside its caller-declared speed-ratio validity
   range.
3. The interval below the explicit stall boundary is non-admissible; the solver
   refuses an intersection there.
4. Every loss resistance is finite, positive, typed, and carries source and
   uncertainty authority. Coefficient uncertainty alone never creates a
   validated-domain card. An attached loss domain has at least one named,
   finite axis and an exact retained source-version identity.
5. A complete operating-point result has exactly one `Certified` interval root
   and no `Possible` root boxes for the nominal declared model.
6. Manufacturer/loss/leakage uncertainty is attached as model-form `Estimate`
   evidence. It is never relabeled as a rigorous physical enclosure.
7. Terminal branch order follows deterministic network traversal, and all
   provenance hashes are order-stable and bind the complete fan curve, source,
   tolerance, fan-bank configuration, recursive network topology, loss data
   — INCLUDING each element's declared `regime_validity`, axis names and
   bounds together, with absence encoded as a state distinct from any declared
   domain — and explicit leakage identity. Equivalently: two `LossElement`
   values that compare unequal never share an operating identity, which
   `loss_elements_that_compare_unequal_never_share_an_identity` checks
   exhaustively over the constructible variants. The `regime_validity` clause
   is called out by name because omitting it silently falsified this invariant
   (bead `frankensim-yq435`): that field is the sole input to
   `regime_audit_cards`, so two networks differing only there reach a reviewer
   as an admitted result and a demoted one while hashing identically.
8. Junction and surface declarations canonicalize index order and reject
   duplicates. Equal junction maxima choose the lowest canonical vertex index.
   That tie-break is delivered by the canonical ascending order plus a strict
   `>` comparison, not by an index guard inside the scan: because every
   candidate is strictly greater than the already-visited incumbent, such a
   guard would be unreachable rather than protective.
9. Every thermal QoI budget has exactly eight terms. Widening a valid upstream
   pressure/flow or efficiency interval cannot shrink the corresponding
   conditional term; changing a requirement or efficiency authority rebinds
   the affected QoI identity even when the nominal scalar is unchanged. The
   temperature-QoI identity binds canonical tetrahedral connectivity and every
   physical vertex coordinate, so geometry-only changes also rebind it.
10. Raw temperature extrema, surface summaries, and margin remain
    `NumericalKind::NoClaim` until an admitted DWR/refinement-to-QoI map exists.
    A conduction residual measured in watts is not converted into kelvin by
    dimensional wishful thinking. Admitting a `DiscretizationReceipt` does not
    relax this: a caller-declared kelvin ladder observation populates the
    budget's discretization term and leaves the numerical certificate alone.
11. As emitted by `extract_thermal_qois`, the margin budget equals the
    junction-maximum budget term-for-term. This is checked with a declared
    discretization receipt present, so the equality is sensitive rather than an
    artefact of both budgets being entirely Unknown. The property is stated at
    extraction, not after `audit_operating_envelope`: the audit issues one
    receipt per declared QoI, so a caller that declares different consumed
    cards for the maximum and the margin can legitimately demote them
    differently.
12. Pressure uncertainty endpoints come from the independently solved
    low-fan/low-resistance and high-fan/high-resistance corners. A flow root is
    never recombined with the opposite resistance bound to manufacture a wider
    but physically unreachable band.
13. Every successfully admitted `AirPath` and every published march or solve
    result keeps its represented conjugate reductions inside the finite `f64`
    domain. Path admission checks the derived positive capacity rate, segment
    conductance, NTU, effectiveness terms, and total wetted area. The march
    checks its temperature difference before products, every segment heat, and
    the running air-heat total. The coupled loop checks solid totals, the
    air/solid difference, each area-weighted residual contribution, the running
    weighted sum, and its final division. Final assembly independently
    rechecks every region difference and both heat-rate totals before it can
    publish a balance. Zero remains legal for signed heat rates, differences,
    and residuals; it refuses where the model mathematically requires a
    strictly positive denominator or exchange factor. The separately
    documented infallible decomposition attachment is a poisoned diagnostic,
    not a published finite-result guarantee.

## Error model

`AirflowError` refuses malformed curves, invalid tolerance or speed domains,
empty network groups, zero/invalid resistances, stall operation, absent curve
intersections, incomplete/ambiguous root searches, unknown branches, and bad
convection-handoff inputs. Loss-card attachment additionally refuses an
unconstrained/unusable validity box or malformed projected identity. Caller
input does not panic. Conjugate path aggregates that must be positive refuse as
`InvalidConjugateInput` with their exact field; signed runtime arithmetic
refuses as `NonFiniteCoupling` with the exact producing stage. A finite row
cannot disappear into a non-finite multi-row total.

## Determinism class

Interpolation and network traversal have fixed order. Square roots use
`fs-math::det`; the numerical root uses deterministic `fs-ivl` subdivision.
Repeated execution is checked against the complete public `OperatingPoint`
artifact, with explicit bit checks on the nominal values and certified root
endpoints. Results are intended to be bit-stable on the same ISA. The solver
has no internal parallel scheduler, so worker-count invariance is not a
separate execution mode. Cross-ISA G5 evidence is not yet retained.

## Cancellation behavior

Curve evaluation and each finite network reduction are bounded scalar work.
The interval search has an explicit 65,536-box ceiling and returns a structured
refusal rather than running without bound. No asupersync cancellation poll is
required for that rung.

`conjugate::solve_conjugate` DOES poll: every outer iteration opens with
`Cx::checkpoint`, so a cancelled exchange stops at an iteration boundary and
returns `AirflowError::Cancelled` carrying the index of the iteration that had
not yet run. Each outer iteration is also a resume boundary:
`ConjugateIteration::reference_temperatures_k` is the complete state that
iteration was solved against, `AirflowError::Cancelled` carries that vector so
a caller can resume from the refusal itself without having retained a history,
and `solve_conjugate_from` restarts from either. The poll happens at iteration
entry AND after the exchange that meets the criterion, so a cancellation
requested during the final solid solve refuses instead of publishing success.
A resumed run restarts iteration numbering, the energy audit, and the history
at zero: it is a fresh exchange seeded from a retained state, not a
continuation of the original record.
Under `Relaxation::Fixed` the resumed tail reproduces the uninterrupted run
bitwise; under `Relaxation::Aitken` the scalar relaxer's own omega history is
not carried across the resume, so the tail is a valid continuation and not a
bitwise replay.

## Unsafe boundary

None. Workspace unsafe-code denial applies.

## Feature flags

None.

## Conformance tests

- G0 fan interpolation, monotonicity refusal, and series/parallel resistance
  composition;
- G0 identical-fan series/parallel affinity identities;
- explicit stall refusal and a sign-changing, unique interval root bracket;
- three declared fan speeds, complete-artifact repeat replay, and
  leakage-resistance sensitivity;
- clause-addressed fan-law authority and a pressure-corner falsifier that
  rejects cross-pairing each solved flow with the opposite resistance;
- typed branch velocity/Reynolds handoff into `fs-convection`;
- G0/G3 plate-fin end-to-end wiring at three fan-solved speed ratios:
  card-backed base-to-fin bondlines, separately attributed fin-flank and
  base-floor Robin heat, an explicit channel-mean bulk-temperature closure,
  a wrong-hydraulic-diameter falsifier, and structured refusal of a
  fan-solved Reynolds point above the laminar-card ceiling. At fixed geometry
  and fluid properties, the developing source-table row makes fin-flank Nu
  and scalar `h` increase with each solved Reynolds point, while the companion
  fully developed limits correctly keep scalar `h` invariant. Evaluation
  provenance binds every changed point;
- semantic-identity separation when uncertainty authority changes without
  changing the nominal operating point.
- G0 sharp-edged-orifice vent lowering: closed-form resistance derivation
  with an independent operation order, allowance coverage of both sourced
  `Cd` endpoints, non-physical-input and empty-name refusals, plateau regime
  card identity/axis/bounds, and a fully vent-lowered enclosure (two parallel
  orifice vents plus an orifice leakage seam) solving to a certified
  operating point with the leakage branch reported.
- actual loss-card projection from an explicitly validated `LossElement`,
  including refusal of unconstrained pseudo-validity and an isolated
  out-of-domain loss Reynolds perturbation that demotes all seven E05.10
  records and rebinds their model-form terms to the exact receipts;
- G0 E05.10 fixture emitting all five QoI families and seven records (the
  uniformity family has mean, spread, and face-mean standard deviation), each
  with a complete eight-term budget and term-by-term provenance rendering;
- deterministic region-order/tie-break equivalence, missing-requirement and
  malformed-region refusals, G3 upstream-envelope widening monotonicity, and
  source-only identity rebinding for fan power and margin.
- G1 analytic surface oracle: for `T = a + b*z` the closed unit cube's
  area-weighted mean is `a + b/2` and the surface-vertex spread is exactly `b`
  by calculus, and the junction maximum is `a + b*max(z)` over the declared
  region. The expectations come from integration, not from re-running the
  production loop. The constant-field companion pins spread and face-mean
  dispersion at exactly zero, and the dispersion is bounded by half the spread.
- G0 margin propagation: with a declared `DiscretizationReceipt` the junction
  maximum's discretization term becomes the declared interval and the margin
  inherits every term. The receipt-absent arm asserting a named `Unknown` is
  the non-vacuity guard: without it the equality would also hold for a
  fabricated budget, because both budgets are otherwise entirely Unknown.
- G0 safety-factor single application: holding the effective limit fixed while
  changing only the declared factor leaves the margin value bit-identical and
  still rebinds the budget identity. A mutant that re-applies the factor fails
  the value assertion; one that drops it from identity fails the identity
  assertion. Factors below one and non-finite factors refuse.
- G3 validity narrowing: holding one operating point fixed and shrinking a
  consumed card's `flow_m3_s` bound around it moves coverage from
  `FullyInDomain` to `FullyOutOfDomain` and demotes every record's model-form
  term. The in-domain arm is the non-vacuity guard. This is the card-narrowing
  direction, distinct from the existing point-moving demotion tests.
- G3 margin offset invariance: shifting the effective limit and every nodal
  temperature by a common offset leaves the margin and the surface spread
  invariant and shifts the mean by exactly that offset.
- G4 upstream-evidence mutation sweep: degrading one upstream slice at a time
  against an admissible baseline. Declared material provenance changes the
  model cards and the `Parameters` reason on every temperature QoI (the margin
  reaches the right reason only by inheriting it, which exercises the
  propagation seam); a receipts-claimed/zero-receipts report refuses; a
  non-converged linear record refuses naming the failing solve, while a
  converged record of the same shape stays admissible, so the guard keys on
  convergence rather than on the record count.
- G1 conjugate air path: the single-segment march against the closed-form
  heated channel `T_w - (T_w - T_in)e^(-NTU)` evaluated through an independent
  `f64::exp` path; the defining identity `h A (T_w - T_ref,eff) == Q` across
  four decades of `h` and three of area; and physical-interval containment up
  to `NTU ~ 25`.
- G3 segment-refinement invariance: splitting one uniform-wall channel into
  `N` equal segments leaves the outlet temperature and total heat unchanged for
  `N` up to 256, EXACTLY, because `e^(-NTU) = (e^(-NTU/N))^N`. The rejected
  arithmetic-mean model is built inside the test and shown to fail the same
  check with a defect that grows with `NTU` (sub-millikelvin at `NTU ~ 0.6`,
  above 10 K at `NTU ~ 20`), so the invariance is demonstrably sensitive rather
  than merely satisfied.
- G1 coupled fixed point against a manufactured solution: over a lumped solid
  the exchange has closed forms for both the outlet air temperature
  (`T_in + P/(m_dot c_p)`, a global energy statement) and the solid temperature
  (`T_in + P/(m_dot c_p (1 - e^(-sum NTU)))`); neither is evaluated by the
  driver.
- convergence battery: monotone residual contraction, Aitken reaching the same
  fixed point, and a typed non-convergence refusal that withholds the reached
  temperatures.
- G4 drills: cancellation at an outer-iteration boundary with the solid side
  proven not to re-run, and checkpoint resume reproducing the uninterrupted
  tail bitwise.
- G0 derived-arithmetic refusals: positive mass-flow and heat-capacity inputs
  whose product underflows, segment conductance/NTU range failures, overflowing
  total wetted area, per-segment and multi-segment air heat, multi-region solid
  heat, the air/solid difference, both area-weighted residual stages, and every
  independently recomputed converged-balance reduction retain exact diagnostic
  attribution.
- flux-balance sensitivity: the per-region audit closes to 2.4e-11 W on a 1.2 W
  fixture AND is shown to open under an injected dropped-face fault that the
  converged temperatures alone cannot reveal.
- G0/G3 conjugate end-to-end (`tests/conjugate_e2e.rs`): a solved fan operating
  point drives a validity-gated `CircularDuctLaminarCwt` card into a real
  `fs-conduction` FEM solve at `Re ~ 854`, `NTU ~ 2.35`, converging in 32 outer
  iterations; the air carries exactly the dissipated 1.2 W, the wall
  temperature rises monotonically downstream with a 14.9 K stream-wise tilt, a
  faster fan lowers peak, tilt, and air rise together, and the exchange replays
  bitwise. A companion test solves the SAME bar under the old one-way declared
  ambient and shows the coupled peak is 8.4 K hotter, so the coupling is not a
  refinement of the previous model but a materially different answer.
- CHT transfer: `restrict o prolongate = identity` on the coarse space, and the
  property `Refine1d` cannot state — prolongating a wall state onto the refined
  path preserves the air-side outlet and heat rate exactly at refinement
  factors 2, 5, and 16.

## No-claim boundaries

- Fan points and tolerances are caller-supplied. Synthetic fixtures prove API
  and algebra behavior, not a manufacturer product's performance.
- AMCA Publication 201-02 (R2011), section 6.5.1, justifies the retained
  speed-scaling relation, not arbitrary operation: the curve's declared
  speed-ratio domain remains binding.
- Repeat replay proves the current scalar implementation is stable on the
  exercised same-ISA path. It does not claim an internal worker-count mode that
  the crate does not implement, or cross-ISA equality.
- The plate-fin fixture is synthetic and one-way coupled. It does not solve
  channel momentum or fluid energy: its mean bulk-air reference temperature
  is the midpoint implied by declared base power, constant density/heat
  capacity, and the solved branch flow. It is neither conjugate CFD nor
  manufacturer or experimental validation. The `conjugate` module supersedes
  that closure for callers who want a two-way exchange; the plate-fin fixture
  is deliberately left as-is so the one-way and coupled answers stay
  separately inspectable.
- The convergence criterion is TWO gates, not one. A reference-temperature
  residual in kelvin cannot bound an interface heat rate in watts without a
  conductance, so the interface-imbalance gate is a separate declared
  obligation and a run that meets the temperature criterion while its worst
  per-region imbalance exceeds it refuses with
  `AirflowError::ConjugateBalanceUnclosed`. The dropped-face fixture proves the
  gap is real: it reaches its temperature fixed point and is caught only by the
  watt gate.
- The watt gate is HYBRID, not a constant. The admitted imbalance is
  `max(balance_tolerance_w, balance_relative_tolerance * scale)` with
  `scale = max_j max(|Q_solid,j|, |Q_air,j|)` from the response being judged,
  and the refusal reports that effective threshold. An absolute constant alone
  is scale-dependent in both directions — a microwatt-scale wiring fault hides
  under a watt-scale constant (measured: the dropped-face fault at 2e-6 W
  produces ~2.5e-7 W of signal), and kilowatt-scale floating-point residue can
  spuriously refuse — and a mW-vs-W unit mistake shrinks signal and gate
  together, which is exactly when the gate is needed most.
  `balance_relative_tolerance` admits `[0, 1)`; `0` restores a purely absolute
  gate. The gate still cannot certify anything smaller than the floor on an
  interface whose true scale is zero.
- Relaxation factors are ADMITTED, not merely finite. `Fixed` requires
  `0 < ω ≤ 2` and `Aitken` requires `0 < ω_init ≤ ω_max`; an inert `ω = 0`
  or beyond-over-relaxation factor refuses at admission as
  `InvalidConjugateInput` instead of stalling for `max_iterations` solid
  solves and blaming convergence. A reference vector the relaxation itself
  overflows refuses as `NonFiniteCoupling` at the "relaxed reference
  temperature" stage, attributed to the relaxation rather than to the next
  solid solve. `ω_max` is otherwise uncapped and is the caller's own risk
  budget.
- `SolidRegionState::mean_reference_temperature_k` is an OPTIONAL wiring
  claim. `Some` is held to the reference vector the driver sent (within
  1e-9 relative, kelvin-denominated) and refuses as
  `ReferenceTemperatureMismatch`; `None` makes no claim and skips the check —
  it is a no-claim boundary, not a wiring proof. `from_robin_flux` always
  claims, because `fs_conduction::RobinFlux` carries the value for free. The
  watt gate cannot see a reference-wiring error smaller than `hA` times its
  tolerance; this check is denominated directly in kelvin.
- `SegmentRefinementTransfer::{restrict, prolongate}` POISON a state whose
  arity disagrees with the bound path: the infallible `fs_ladder::Transfer`
  contract cannot refuse, so the mismatch returns all-NaN at the correct
  output arity rather than silently truncating the downstream-most entries or
  restricting a partial block with full-block weights. Downstream response
  gates refuse non-finite walls.
- "Exact" in this crate's conjugate documentation means exact in REAL
  arithmetic. Measured floating-point bands on the current fixtures: segment
  refinement holds the outlet to 1e-10 K and the heat rate to 1e-8 W up to
  N = 256; the transfer round trip lands within 1e-9 K rather than bit-exact,
  because the restriction evaluates a ratio of geometric sums; the coupled
  per-region balance closes to 2.4e-11 W on a 1.2 W fixture.
- `SegmentRefinementTransfer::restrict` is the DOWNSTREAM-WEIGHTED mean
  `W = (1-a) * sum_i a^(f-1-i) W_i / (1 - a^f)` with `a = e^(-NTU/f)`, not a
  plain mean. A plain mean is outlet-preserving only on the uniform blocks
  `prolongate` emits, so a round-trip test cannot distinguish them; on a
  nonuniform fine state the plain mean is wrong by 15.65 K at NTU 2 over four
  sub-segments. The nonuniform check in `tests/conjugate.rs` is the one that
  discriminates.
- The conjugate exchange is NOT CFD. Air is a stream-wise chain of well-mixed
  1-D segments: no lateral mixing, recirculation, buoyancy, heating-driven flow
  redistribution, or momentum coupling back to the operating point. The
  `RANS` and `LES` rungs remain declarations.
- `h` is FROZEN across the outer loop. The coupling variable is the reference
  temperature vector alone; air properties are not re-evaluated at the drifting
  film or bulk temperature, so a temperature-dependent `h` is outside the
  model. A caller who wants that must re-run the driver with a new `AirPath`.
- Relaxation is SCALAR. `fs_couple::AitkenRelaxation` is a scalar delta-squared
  relaxer; the driver projects the vector residual onto one area-weighted
  scalar and applies a single omega to the whole vector. This is not a vector
  interface accelerator. On this seam the composite gain is
  `(1 - eps/NTU)` times the solid's gain, both strictly below one, so plain
  staggering already contracts — the retained Aitken path is wiring, and no
  claim is made that relaxation rescues a divergent case, because the battery
  contains none.
- The per-region flux balance is a WIRING falsifier, not a conservation proof.
  Once `T_ref,eff` is defined so that `h A (T_w - T_ref) == Q`, the per-region
  balance is an algebraic identity at the fixed point for a uniform `h`; it
  catches a dropped face, a mis-bound region, or a wrong mass flow, and it
  catches nothing about physics. The independent cross-check is
  `decomposition_residual_w`, which compares the declared regions' sum against
  `fs-conduction`'s own whole-domain `robin_out_w` loop and so detects a Robin
  face owned by no region or counted twice. That is also wiring. The physical
  content lives in the closed-form and refinement-invariance checks above.
- `ConjugateSolution::with_decomposition_cross_check` is an INFALLIBLE
  compatibility surface. It cannot return a structured refusal without a
  public API change. A non-finite whole-domain Robin total or subtraction
  poisons `decomposition_residual_w` with `NaN`, which withholds a finite
  diagnostic but does not identify which operand was invalid. Callers needing
  exact refusal attribution must validate both operands and ensure their
  subtraction is representable before attaching the total; this method alone
  is not an admission gate.
- `AirSegment::ntu(capacity)` is also an INFALLIBLE compatibility helper. It
  performs unchecked raw arithmetic for an arbitrary caller-supplied capacity
  and can return zero or a non-finite value. The checked `AirPath` admission
  and march paths validate capacity, conductance, NTU, and effectiveness before
  using them; callers must not treat the standalone helper as an admission
  receipt.
- The seam is typed, not ledgered. `conjugate::SEAM_PORT_KIND` declares the
  `fs-couple` port kind and its effort dimension is checked against
  `Temperature::DIMS`, and `fs_couple::EnergyAudit` records every exchange's
  imbalance with NaN poisoning. The window-balance ledger path
  (`BoundaryTemperatureReference`, `WindowEvidenceRef`) is NOT wired, so no
  ledgered entropy accounting is claimed.
- `SegmentRefinementTransfer`'s fine side is a REFINED CORRELATION STATE, not a
  RANS field. It defines the correlation rung's own refinement semantics and
  the state shape a RANS rung would have to accept; it is not evidence that a
  RANS rung exists. The `RANS` and `LES` relative-cost hints in the CHT ladder
  remain unmeasured declarations, because those rungs do not exist to measure.
- The conjugate fixtures carry no maturity registration, retained corpus
  receipt, machine fingerprint, or experimental comparison. No capability
  level may be inferred from them.
- The fin-flank row is the narrow Shah-London Chapter VII, Table 52 slice for
  smooth simultaneously developing rectangular flow at `alpha*=0.5` and
  `Pr=0.72`. Its lower-`Gz` bridge is a declared engineering interpolation,
  not a published curve. It does not model a commercial plate-fin array,
  shroud/bypass flow, conjugate heat transfer, or interrupted-fin effects.
- The companion rectangular CWT/CHF rows are fully developed analytic limits.
  Reynolds is a validity axis but not a formula variable for those rows, so
  their `Nu` and scalar `h` stay invariant when only fan speed changes.
  Evaluation provenance still binds the changed Reynolds point, and a
  sufficiently high fan-solved point refuses rather than extrapolating.
- Summing named Robin heat rates back to the whole-domain Robin total is only
  a boundary-wiring falsifier; both reductions use the same face quadrature.
  The wrong-`D_h` case likewise proves the geometry handoff is load-bearing,
  not that any hydraulic diameter produces validated hardware physics.
- Piecewise-linear interpolation and quadratic loss coefficients do not model
  swirl, recirculation, acoustic interaction, thermal buoyancy, compressibility,
  fouling, or installation system effects.
- A `LossElement` with no attached validity remains an estimate with no regime
  card. An attached domain is only the source owner's stated operating box; it
  does not validate the quadratic loss law or turn its uncertainty into a
  rigorous enclosure.
- The sharp-edged-orifice lowering models a clean thin sharp-edged opening
  only. It is not a discharge model for louvered, screened, filtered, or
  ducted vents, claims no compressibility or installation effects, and its
  declared allowance is an engineering spread over sourced `Cd` values, not a
  published confidence interval. Density is a caller declaration; deriving it
  from an operating envelope is the lowering stage's obligation, and no
  retained orifice measurement corpus exists.
- The nominal root bracket certifies only the declared mathematical model.
  Tolerance-propagated flow, pressure, branch splits, and Reynolds values remain
  `Estimate`, not validated hardware envelopes.
- Parallel fans are identical and equally loaded. Unequal curves, unstable
  parallel operation, active control, fan-fan interference, and transient
  startup remain outside this slice.
- No retained manufacturer table, wind-tunnel corpus, CFD comparison, or
  experimental enclosure validation exists; there is no L4 or product claim.
- The thermal QoI consumer does not close E05.10's external validation, DWR,
  mesh-refinement, sensor, or naked-scalar lint obligations. It emits the rich
  budget and now enforces final validity intersection when the orchestrator
  supplies the complete consumed-card registry, card-use map, and operating
  envelope. Completeness/authenticity of those supplied authorities remains an
  orchestration and package/checker responsibility; the broader E05.10 bead
  remains open.
- The conditional pressure/flow and efficiency intervals are only as sound as
  their caller-declared source envelopes and the stated quadratic model. The
  always-explicit unknown model-form term prevents their interpretation as a
  whole-product uncertainty bound.
- Fan power means input power under the declared total efficiency,
  `Delta p * Q / eta_total`. Motor/controller transients, reactive power,
  acoustic power, and installation effects are outside this slice.


### Porous Darcy and mixed passive pressure graphs

`graph::LossGraph::new_mixed` admits typed linear Darcy loss, existing quadratic
loss, and an explicitly authored sum `dp = A Q + B Q |Q|`. The series/parallel
quadratic API keeps its original meaning. `graph::porous::LinearLossElement`
can bind the shared material barrier and coherent gas viscosity; layer resistance
is compiled once and graph iterations use numeric coefficients. A sealed path
is omitted, and an inertial/Forchheimer coefficient is never inferred from k.

`PassiveBranch::from_darcy_forchheimer` binds both coefficients to the same
original barrier, area and coherent gas snapshot. The explicit superficial-
velocity calibration supplies `B = rho sum(beta L)/area²`, while the original
intrinsic tensor reduction supplies `A = mu sum(L/k)/area`. Graph admission
checks the retained pair; dropping the calibrated quadratic term from the
retained linear binding or changing its value refuses. Exact zero inertia
recovers a linear branch. The existing fan, conservation and sensitivity owners
consume these coefficients without a second graph solve or material database.

Material edges share one coherent frozen gas state. Its absolute pressure is
the common reference at graph gauge pressure zero. Actual endpoint absolute
pressures and their excursions from that reference must remain inside source
support and the caller's constant-density allowance; superficial speed is
checked at the accepted solution. The compressible finite-volume material
transport owner remains the path for substantial density change or gas-specific
slip. This graph does not resolve membrane diffusion, rarefied pores or storage.

Mixed fan intersections solve the actual pressure graph at each trial instead
of applying quadratic homogeneity. Numerical brackets may leave a material's
operating domain, but final acceptance cannot. Tangents and adjoints use
`dQ/ddp = 1/(A + 2 B |Q|)`; a branch's resistance-scale derivative scales both
coefficients together. Linear stagnant edges have their finite physical
derivative; pure quadratic stagnation retains its existing refusal. Analytic
and manufactured-network regressions cover reversal, extreme coefficients,
series/parallel/loop balances, gas substitution, source admission, fan operating
points and sensitivity agreement. These remain nominal model estimates, with
source uncertainty retained rather than promoted to a numerical certificate.

### Gravity in material pressure networks

`LossGraph::with_gravity` binds caller-supplied Cartesian node positions and
signed gravitational acceleration to one actual `FluidTransport` snapshot.
There is no Earth-gravity default. Each edge obeys
`p_from - p_to + rho g·(x_to-x_from) = A Q + B Q |Q|`; gravity is a conservative
pressure supply, never a resistance or a heat source. The solver works in
pressure relative to the uniform reference-density hydrostatic field, so a
quiescent column or loop balances without numerical relaxation. Source queries
and head construction happen once at binding; hot iterations use the existing
cached resistances and cached head coefficients. Node and edge traversals keep
the original cancellation and deterministic-order behavior.

`with_boussinesq_gravity` additionally binds one source density per edge and an
explicit bound on `abs(rho-rho0)/rho0` in `[0,1)`. Only density departures drive
the reduced pressure equations. This enables frozen hot-leg/cold-leg natural
circulation while retaining ONE reference density for volume/mass conservation.
The original source snapshots remain inspectable. Drag coefficients are not
silently reinterpreted; updating temperature, composition or material requires
rebinding the affected source/loss descriptions. Bound Darcy gas density must
agree with the reference density. This approximation does not provide a
variable-density mass solve, self-consistent compressible atmosphere, evolving
thermal circulation or experimental validation.

Every accepted branch retains static pressure drop, full gravitational head and
the actual dissipative drive separately. `dissipative_pressure_drop * flow` is
the passive loss; static pressure alone may rise downstream. Fan intersections
use actual pressure-graph trials when gravity is present, because quadratic
homogeneity no longer applies. Existing tangents/adjoints hold gravity,
positions and density fixed, differentiating the dissipative drive rather than
mistaking conservative head for a resistance coefficient.

The existing `TransportNetwork` consumes gravity-driven flows, including
downward paths with increasing static pressure, using the same reference rho/Cp.
The original `TransportNetwork::new` retains directed-acyclic admission.
`new_circulating` additionally transports heat around a closed native hydraulic
loop with at most 64 nodes. All pressure boundaries are explicit pressure
datums; incoming and outgoing flow must balance there as well as at interior
nodes. Each active component needs a numerically resolved wall-heat anchor.
An entirely adiabatic closed component has no unique steady temperature and
refuses. Stagnant adiabatic nodes retain their unknown temperatures.

The same exponential `AirPath` laws and incoming-capacity mixing weights form
a simultaneous node-temperature system. The existing bounded elimination and
backward-residual check solve that system; an independent physical branch
remarch checks the sum of absolute node heat residuals and the original
whole-network heat budget. The supplied positive temperature reference only
sets numerical coordinates and the initial solid-iteration guess. It supplies
no fresh air or thermal reservoir. Uniform-wall equilibrium is evaluated in
wall-centered differences to retain the constant-temperature limit at weak hA.

Fixed-flow wall-temperature/hA tangents and adjoints close through the same
node matrix and its transpose; conductance scaling preserves this mode.
Imposed external return-air links and individual unconstrained branch-flow
derivatives refuse on a closed circulation. Density, heat capacity, branch
flows and hA remain frozen. There is no temperature-driven hydraulic feedback,
fluid storage, conditioning enclosure or experimental validation claim.
Neither gravity work nor hydraulic dissipation is automatically deposited as heat.
Focused G1/G3/G4 regression sources compare actual buoyant hydraulic loops with
independent exponential heat balances, check reversal, numerical-reference
invariance, weak-anchor uniform equilibrium, derivatives, refusal and cancellation.
G0/G1/G3/G5 focused tests cover arbitrary-axis hydrostatics, zero/reversed g,
material-density substitution, closed-loop density-driven circulation and power
closure, fan operating points, fixed-head derivatives, heat-transport coherence,
invalid inputs and deterministic repetition. These are analytic-model tests,
not a measured performance or physical-validation claim.


### Permeability with mixed finite-gas heat laws

`finite_gas::network::exchange::transport_paths_step_candidate` and
`advance_with_transport_paths` pass the complete material `GasTransportPaths`
through the original `MaterialExchangeState` owner. Darcy leakage and selective
membrane changes therefore reach conduction, forced/natural convection, reaction,
finite storage and wall motion through the actual transported gas. Each stage
runs once, carries partial enthalpy once, and publishes every endpoint together.
Final transport driving-force checks enclose the thermal/mechanical solve; the
constitutive source domains, prescribed convection geometry and first-order
splitting limits remain explicit. Existing two-mechanism entry points delegate
to this same owner.


## Stratified wall-opening exchange for compartment networks

`stratified::wall_vent_flow` integrates the original quadratic orifice law
between two hydrostatic, two-layer columns. It partitions each rectangular
opening at both layer interfaces and every neutral-pressure plane, producing
up to six directional slabs. Each slab retains its donor and receiver layer,
uses the actual donor density, and analytically integrates the square root of
the linear pressure difference. Small counterflow intervals are evaluated in
local coordinates rather than snapped away when world coordinates round alike.

The caller supplies a sourced discharge coefficient and the original resistance
uncertainty semantics. Explicit Cd lowers once to `LossElement`; the convenience
sharp-edged path calls the existing Cd=0.61 owner unchanged. CFAST's Cd=0.7 is
not silently substituted. Source identities and validity declarations persist;
this leaf does not propagate uncertainty or certify an empirical regime.

Returned slabs preserve simultaneous counterflow even if the net mass rate is
zero. Layer mass rates are equal and opposite across the opening; optional
enthalpy advection uses each slab's donor enthalpy and the same mass transfer.
The primitive does not mutate species or thermal inventories. Zero area seals
exactly, zero gravity is supported, and malformed/unrepresentable columns or
cancellation refuse without a partial result. Successful evaluation allocates
no vectors or per-slab state.

Ten focused native tests pass: sealed/equal columns; original orifice equivalence;
donor-density reversal; analytic neutral planes; all six slabs; opening
refinement; zero-net-mass heat exchange; reference-height invariance; tiny
counterflow; and invalid data/cancellation. Tests compile the actual new source
and exact original loss/Cd/error owners with the original pure dependencies in
an isolated leaf build, not the entire airflow/conduction/runtime graph.
Reference: [NIST TN1889v1, section4.1, equations4.1–4.4](https://doi.org/10.6028/NIST.TN.1889v1).
This remains a low-Mach wall-opening closure; ceiling exchange, jet entrainment,
compressible choking and experimental validation require separate owners.


## Conservative compartment-fire networks and thermal walls

The optional `compartment-fire` feature composes the existing source-owned
`fs-gas::multicomponent::compartment` rooms, stratified wall-opening owner,
`fs-la::prepared_lu`, `fs-thermal::enthalpy`, `fs-convection::fire`, and the
original `fs-conduction` surface-emissivity/view-factor owners.
The bounded network admits up to 32 source-matched rooms and 512 wall, ceiling,
or mechanical connections. Each directional vent slab carries its actual donor's
species and original partial enthalpies, with identical opposite room entries.
Ambient exchange and explicitly supplied fan work enter the external account.
Openings may be changed or closed exactly between accepted steps.

A bounded damped Newton solve treats opening pressure implicitly while donor
temperature, composition and interface geometry are frozen per attempted step.
All trial rooms are reconstructed from accepted state. Whole-network duration
retries and cancellation publish no partial rooms. Duration-dependent source
providers recompute wall candidates once per actual attempted interval; only
the accepted payload is returned. A subsequent split-source provider operates
on transported inventories before publication, preserves room identity and
clocks, and declares its actual additional sources. The report distinguishes
post-source pressure changes from the transport pressure-equation residual and
reports net physical external energy, each species source, and balance residuals.
This remains a first-order estimated operator composition, not a globally
second-order or fully temperature-implicit scheme or discrete entropy proof.

`ZoneNetwork::step_with_pressure_sources` also evaluates a finite source inside
every original Newton, Jacobian and line-search pressure trial. Each evaluation
first reconstructs the accepted inventories, applies transport, then supplies
those actual transported inventories to the source callback. The pressure
equation closes against the source-updated rooms, so chemical heating affects
opening flow in the same accepted step. It does not burn oxygen that has already
left through a vent. The original Newton/LU, identity/clock checks, balances and
transactional publication are shared with the existing entry points.

Preparation data are immutable during each attempted interval; trial diagnostics
and candidate external states are returned separately, with no cloning
requirement. Only the converged trial's result is returned. A failure or
cancellation after a chemical update leaves all accepted rooms unchanged, and
duration retries rebuild preparation. Its pressure residual includes the source;
the post-pressure-solve split increment is zero because there is no subsequent
source stage. Donor temperature, composition and interface geometry in opening
fluxes remain frozen at the start of the step. This pressure closure does not
establish full thermochemical implicitness, second-order accuracy, equilibrium
chemistry or an entropy theorem. Eight focused regressions cover exact sealed
oxygen capacity, post-outflow availability, an independent reconstruction of
opening transfer at the accepted chemical pressure, multi-room redistribution,
retry/cancellation atomicity, source admission and zero-source equivalence.

The wall adapter retains the caller's original borrowed ThermalLaw, original
DepositionSimulation histories and BackwardEulerStep residual/acceptance gates.
Uniform through-thickness finite volumes use original series conductances and
both surface half-cell resistances. The actual material masses and geometry
remain fixed when layer exposure areas change. Signed gas-contact energy is
returned identically opposite to the wall exchange; outer convection and
explicit incident radiant power have separate accounts. Incident fire radiation
is external wall input, not an additional gas debit. Gas emissivity, soot
absorption, thermal properties and inverse calorics are never inferred.
Numerical vectors/LU are reused; original thermal trial maps still allocate
within explicit cell, Newton and backtracking budgets.

`Wall1d::candidate_with_exchange` accepts a separate signed internal radiant
power at the inner face. Negative power represents a radiating wall losing
energy. This power enters the original half-cell surface balance and thermal
storage residual, and is reported as `coupled_radiant_energy`; it never becomes
`external_radiant_energy`. The caller applies the matching enclosure gas and
other-wall powers exactly once. The old `candidate` entry point delegates with
zero coupled radiation. Holding radiation fixed over the interval does not
establish unconditional stability or a fully implicit radiation/temperature
solve; coupled iteration or timestep refinement remains required.

The prescribed-fire helper computes lower-layer density and mass-specific Cp
from the original current species/caloric inventory, then reuses Heskestad
entrainment. One identical donor vector moves from lower to upper; convective
heat enters the upper layer and radiant power is returned separately. An absent
upper layer can be created by real entrainment. No fuel/products/soot are
invented, and this prescribed-HRR owner does not consume oxygen. Pyrolysis and
field-resolved flame behavior remain outside this prescribed-HRR
helper. The separate inventory-limited burning
adapter consumes actual reactants and accounts for original formation energy.

`PrescribedFire::evaluate_with_penetration` optionally applies
`PlumePenetration::NeutralBuoyancy`. NIST TN 1889v1 §3.3.3 equations 3.18–3.19
give the equal-composition constant-Cp limit `Qc/[Cp(Tu-Tl)]`. This extension
uses the current lower/upper molecular weights to find the donor temperature
at the upper-layer density, `Tn = Tu R_upper/R_lower`, and integrates the
original donor-species NASA enthalpy differences to that temperature. Accepted
entrainment is the smaller of the original Heskestad estimate and `Qc/Δh`.
An absent upper layer or an already buoyant donor imposes no barrier. A neutral
state outside the original caloric domain refuses rather than extrapolates.
The original `evaluate` retains its unconfined choice. A source with zero
convective power produces no gas heat or entrainment and reports its radiant
power without invoking the flame correlation.

This is a steady reduced penetration limit, not resolved momentum overshoot,
detrainment, or a stalled-plume calculation. Supplied convective heat still
enters the upper layer once. Nominal and accepted entrainment, neutral
temperature, required enthalpy rise, and whether the cap was active are
returned separately. Exact species transfer, common-pressure energy closure,
the published constant-Cp limit, variable-Cp/unequal-composition cases,
caloric-reference and species-order invariance, and domain/cancellation
refusals are tested. No improvement in experimental bias is inferred from
this additional model option.

Focused native verification passes 30 original-source airflow tests: ten
stratified vents, seven network/source-composition cases, nine thermal-wall
cases and four plume cases. These include simultaneous counterflow, independent
species/enthalpy balance, pressure equilibration, forced complete-step retries,
matched gas/wall clocks, joint gas-wall energy, latent heat storage, surface
resistance, implicit conduction beyond the explicit limit, and cancellation.

### Horizontal ceiling/floor openings

`compartment::ceiling` implements the Cooper/CFAST shallow horizontal-opening
correlation from NIST TN 1889v1 section 4.2. It retains both nonnegative
directional flows. Denser gas above lighter gas produces exchange at zero
pressure difference; imposed pressure progressively suppresses that exchange
through the published flooding threshold. Each mass flow uses its actual donor
density. The original discharge/loss owner retains the published coefficient
0.68 and the caller's explicitly declared engineering resistance allowance.
Round and square apertures use the separate characteristic lengths in the
pinned CFAST implementation; other shapes are not guessed.

`ZoneNetwork::new_with_ceiling_vents` places these flows inside the same implicit
pressure solve as wall openings and uses the same species/partial-enthalpy
transport account. A room below the aperture connects at its ceiling; a room
above connects at its floor. Ambient may occupy either endpoint. Admission
allows only coordinate roundoff at the shared boundary, and refuses datums
that cannot resolve room height. The network chooses the source-owned adjacent
present layer so an exactly absent layer is not selected through rounding.
Closure and geometry updates are atomic. This remains an estimated shallow,
low-Mach aperture law: thick shafts, choking, plug-holing and vent-induced
plume entrainment are not represented. Focused law and multi-room tests cover
counterflow, pressure suppression, original-source conservation and cancellation.

### Participating gas-layer radiation

`compartment::radiation::ParticipatingEnclosure` combines original admitted
view factors and card-backed surface emissivities with explicitly supplied
reciprocal optical paths through the two gas layers. The path weights must
close the original area-weighted view factors; they are never renormalized to
hide missing power. Ordered isothermal segments use the nonscattering gray
Beer--Lambert formal transfer solution. Expanding that solution gives symmetric
pair-exchange coefficients; the original prepared LU solves all diffuse wall
reflections. Pair powers enter their two thermal owners with opposite signs.

`evaluate_room` uses actual layer temperatures and refuses nonzero optical path
length through an absent layer. `RadiationReport::add_gas_heat` adds only the
gas powers to original `LayerRates`; corresponding signed surface powers enter
`Wall1d::candidate_with_exchange`. A `step_with_sources` provider can construct
the radiation and wall candidates for each attempted duration, add gas
radiation and convective contact once, and return the wall candidates only
when the network step is accepted. External fire radiation remains separately
booked. Reused dense/chord buffers bound the solve; source-bearing result
objects allocate in proportion to the number of boundary patches.

The caller supplies absorption coefficients, their source/domain assessment,
and optical geometry, and must regenerate paths when the interface moves.
This is an estimated gray optical quadrature, not a gas spectroscopy database,
RADCAL, soot-yield model, scattering solver or geometry generator. Numerical
closure and algebraic residuals do not certify optical quadrature or physical
accuracy. Tests include transparent diffuse plates, black gray-medium slabs,
hemispherical angular refinement against the independent exponential-integral
solution, ordered two-layer transfer, opaque shielding, exact equilibrium,
self-view geometry and actual gas-plus-wall storage conservation.

`radiation::soot::RayleighSoot` supplies a constant-index, dilute Rayleigh-particle
absorption closure for caller-sourced refractive indices, particle density,
loading and spectral-temperature domain. It converts actual particulate mass
concentration to volume fraction; soot is never added to ideal-gas species,
moles or pressure. The spectral law is the pinned RADCAL small-particle law.
Its analytic Planck mean uses `4 ζ(5)/ζ(4)` and the exact defining SI constants,
and is independently checked by direct wavelength quadrature. Zero soot or a
zero absorptive index is exactly transparent. Coefficient evaluations allocate
nothing; converting a report into `LayerAbsorption` explicitly records its
particle population, optical source and blackbody weighting temperature.

This weighting temperature is a radiation-spectrum input: the local particle
temperature gives an emission mean and an illuminating blackbody temperature
gives that source's incident mean. One gray coefficient does not exactly
represent both different spectra or spectral filtering along a finite path.
Small particles, dilution, negligible scattering and constant optical indices
remain source assumptions. Using the existing one-temperature gas heat account
also assumes particle/gas thermal equilibrium. No soot yield, transport law,
molecular gas spectroscopy or universal optical constants are inferred. Tests
exercise the original enclosure's signed surface/gas exchange at four soot
loadings as well as the independent spectral integral.

### Wavelength-dependent participating radiation

`radiation::spectral` extends the same enclosure with sourced, contiguous
`SpectralBand` absorption and optional original `RayleighSoot` populations.
Every surface and gas layer emits at its own actual temperature. Soot absorption
is evaluated at each wavelength, retaining spectral filtering instead of
substituting a single Planck mean. Original reciprocal chords, gray surface
emissivity cards, diffuse reflection equations and prepared LU still own the
transfer. The caller explicitly names the assumption extending those total
emissivities over wavelength; no molecular coefficients are inferred.

`SpectralOptics::new_with_surface_bands` alternatively accepts independent
`BandSurfaceEmissivity` sources on the same gas/surface wavelength partition.
These are spectral properties with explicit temperature domains, not relabeled
total material emissivities. The opaque diffuse law uses the same band value
for emission and absorption and its complement for reflection. Every queried
wall temperature is checked in every band. Supplying gray cards together with
these sources refuses ambiguous ownership; the retained optical report contains
the applied spectral sources and an empty gray-card vector. The shared finite
wall/reactive-surface adapter selects these sources without a dummy gray query.
Constant band optics still use one LU factorization per band. G1 tests compare
selective plate exchange against independent Planck integration with unlike
source temperatures; G3 covers exact gray recovery with gas/soot, equilibrium,
and conservative finite-wall feedback. The source must justify the constant
band approximation. Zero band values admit perfect diffuse reflectors: their
material temperature contributes no emission or spectral-coverage restriction,
and reflected irradiation is retained without thermal absorption. Before
factorization, every connected optical component must reach an absorbing
surface or participating gas. An unanchored lossless component refuses even
when another disconnected cavity emits, rather than relying on roundoff to
make its singular system appear solvable. Positive-emissivity gray paths skip
this additional component check. Transmission and directional surface optics
remain outside this opaque diffuse lane.

Positive 8/16-point Gauss--Legendre integration in log wavelength splits at
every supplied band boundary, with initial panel widths at most one logarithmic
unit. Controls bound absolute heat/radiosity error estimates, evaluation count,
and omitted blackbody fraction. Constant band optics reuse one factorization
per band. Reused buffers avoid allocation per wavelength or subdivision.
Reports retain individual output estimates, source identities, wavelength
coverage, omitted blackbody emission, algebraic and energy residuals, and
actual work counts. `evaluate_room_spectral` enforces original layer presence;
`SpectralRadiationReport::add_gas_heat` adds the integrated gas powers once,
while matching surface powers remain available to the original wall owners.

The quadrature differences are estimates, not certified bounds. Omitted
blackbody emission is not a bound on omitted net heat. The source must justify
constant coefficients within each band, optical geometry, negligible scattering
and particle/gas thermal equilibrium. This is not a molecular spectroscopy
database, soot-yield model or RADCAL replacement. Planck normalization retains
the original conduction owner's Stefan--Boltzmann convention.

Six new test groups check independent Planck/tail quadrature, finite-window
gray recovery and factorization reuse, selective absorption bands, independent
spectral Beer--Lambert traversal, equilibrium/adjacent-temperature exchange,
and bounded-work/cancellation/reuse behavior. The existing
`soot_layer_radiation` example accepts an optional spectral suffix while
preserving its original invocation. Its authored comparison and actual native
results are in `data/reference/compartment_fire/spectral_soot_example.json`;
both models conserve energy but predict different layer/wall heat distributions.


### Implicit participating radiation and finite wall storage

`radiation::walls` couples the original gray or wavelength-dependent enclosure
to original `Wall1d::candidate_with_exchange` candidates. Its nonlinear unknowns
are signed surface radiation powers. Each trial advances the finite walls with
those powers and evaluates the existing optical owner at their resulting face
temperatures. The emissivity callback must return the original material receipt
for the actual queried temperature. A previous power vector is only a numerical
initial guess, never a separately applied source.

The prepared workspace reuses original wall, radiation, spectral and LU owners.
A Jacobian column perturbs one wall and the shared optical solve, avoiding a
complete wall-bank advance for every column. Explicit controls bound nonlinear
iterations, wall/optical work, face-temperature and power mismatches, and energy
residuals. Workspaces can be reused after cancellation or rejected candidates;
the complete solve is not claimed allocation-free or unconditionally convergent.

An accepted candidate uses the identical radiation power vector for fresh wall
candidates and the original gas account. Radiation is evaluated again at the
returned wall faces to check the nonlinear residual; the report retains those
differences instead of redistributing them into gas energy. Its `add_gas_heat`
adds radiation plus the opposite wall convection exactly once. A network source
provider installs the returned walls only alongside an accepted gas candidate.
Area, layer presence, wall/room clocks and convection-temperature ownership are
checked. Typed errors preserve cancellation through `NetworkError`.

Gas temperatures, optical coefficients and paths stay frozen within an interval.
Callers regenerate geometry and optical data from accepted gas states. This is
implicit wall/radiation feedback with first-order gas coupling, not a fully
implicit reacting-gas solve, a global second-order method or an entropy theorem.
All existing optical/source assumptions and omitted-spectrum qualifications
remain in force.

Five test groups cover an independent stiff black-plate solution, equilibrium,
selective spectral bands, actual NASA gas/wall storage with forced network retry,
and bounded-work, cancellation and source-accounting behavior. The registered
`radiating_compartment` example evolves two hot gas layers and two finite plates.
Its three actual gray/spectral/refinement runs and conservative storage accounts
are retained in `data/reference/compartment_fire/radiating_compartment_example.json`.
The authored constant-band optics exercise the finite-window gray limit; they
are not a molecular absorption dataset or an experimental validation case.

`radiation::walls::closed_room::candidate_closed_room_walls` adds endpoint
feedback for fixed-inventory, unvented gas and finite walls. Gas heat powers
are the one or two active unknowns. Each trial uses the original `TwoZone`
advance, including its pressure/interface work, and a caller-supplied original
wall candidate evaluated at that gas endpoint. The caller rebuilds opacity,
paths and convection from that supplied state. Bounded Newton/backtracking
reuses a small prepared LU; it never changes the requested timestep.

At acceptance, the gas receives exactly the heat in the returned wall report.
A fresh endpoint candidate checks heat and gas/wall-face temperature mismatch;
combined storage is checked against original external wall heat without any
residual redistribution. Candidate clocks and absent-layer heat ownership are
checked. Failure/cancellation leaves both accepted inputs unchanged. This is a
first-order closed-room heat integrator, not a ventilated/reacting-gas or
unconditional convergence claim. The `radiating_compartment implicit` mode
runs it with actual evolving NASA gas and wall storage. G1 tests compare stiff
radiation/convection with a separately reduced scalar energy equation; G3 tests
two-layer pressure/geometry feedback and time refinement; G4 tests cancellation,
work budgets and incompatible clocks.

`radiation::walls::closed_room::candidate_split_room_walls` applies that same
endpoint thermal operator to an unpublished transported gas candidate through
the original `TwoZone::apply_split`, preserving its clock. Inside
`ZoneNetwork::step_with_pressure_sources`, the returned gas heat enters the
additional source account once and the matching walls are installed only with
the accepted pressure candidate. Radiative cooling/heating therefore changes
vent pressure and flow within that solve. Ambient species and physical stream
enthalpy remain external; internal gas/wall heat cancels in the combined
storage balance. Every trial recomputes opacity, paths and convection from its
actual gas endpoint. Donor vent temperature, composition and interface retain
the existing frozen-start approximation. This is first-order transport/thermal
composition; it does not assert simultaneous chemistry or radiative transfer
through an unrepresented opening. `radiating_compartment ventilated` exercises
the original vent-pressure, species, caloric and finite-wall owners together.
Focused tests independently reconstruct vent flow at the radiatively cooled
accepted pressure, check external physical enthalpy and combined storage, and
check retry/cancellation without publishing a rejected gas or wall candidate.

### Conservative solid-to-gas volatile transfer

`compartment::pyrolysis::VolatileMixture` binds an explicit volatile mass
composition to the complete original `MixtureGas` species bank, caloric domain
and energy references. Its `VolatileEnthalpy` implementation uses the original
frozen NASA-mixture evaluator, including formation enthalpies. A reactive solid
must use this same model for its outgoing stream and must supply virgin/char
calorics with compatible physical formation references. An independently fitted
effective heat-of-pyrolysis chart is not silently interchangeable with that
physical gas enthalpy. No fuel composition is inferred from a material name.

`release_rates` preserves both time-integrated volatile mass and enthalpy from
the original `fs-thermal` pyrolysis owner, including releases produced over
different temperatures. It maps the mass onto the original canonical species
axis and subtracts those species' original internal-energy references to form
`LayerRates`. It does not substitute a final solid temperature or add another
heat of pyrolysis or combustion. The complete receiver gas must match the bound
source. The carried mean enthalpy must remain inside its original caloric
domain; actual rate-integrated mass and reconstructed physical enthalpy must
satisfy explicit numerical gates without repairing the source integrals.

The existing `ZoneNetwork::step_with_sources` provider prepares a solid
candidate for each attempted interval and adds these rates to the receiving
layer once. They already enter the original network pressure equation. The
caller installs the returned solid candidate only after network acceptance;
rejected or cancelled attempts publish neither owner, and shorter retries
recompute the release from the accepted solid. Subsequent gas chemistry uses
the existing source-inclusive pressure-stage API. A raw thermal release record
cannot authenticate its material model, so correct per-cell model binding and
routing remain explicit caller obligations.

Seven integration tests cover original NASA enthalpy and reference conversion,
variable-temperature integrated release, source/domain refusals, cancellation,
joint solid/gas mass and physical-energy closure, and a forced network retry
with source-inclusive pressure. This supplies conservative transfer, not a
volatile-composition predictor, soot model, combustion closure or calibrated
material dataset; those models remain separately owned.

The `pyrolyzing_compartment_fire` example now runs the complete path from an
externally heated spatial solid, through original NASA volatile enthalpy, into
oxygen-limited reacting ventilation. Solid conduction/decomposition determines
the released fuel; no fuel-feed or heat-release history is prescribed. Actual
gas-to-face convection is an internal equal/opposite exchange, and the external
heater/300 K radiative environment has a separate energy account. Every pressure
trial applies original inventory-limited burning after transport. Accepted
solid/gas candidates advance together; whole-step retries rebuild release at
their actual duration. Combined mass, canonical species, C/H/N/O elements,
physical energy, pressure and clocks are checked during execution.

Its four authored small-box cases are closed, restricted/open doorway, and
restricted doorway with pressure-dependent exhaust fan. Exhaust shaft heat
enters downstream ambient and is not counted again inside the room. The solid
is explicitly manufactured to release methane with compatible formation
energy; it is not a PMMA or wood dataset. The case uses one active gas layer,
instantaneous complete oxidation and a simple external optical boundary, not
an ignition, finite-rate flame or smoke-feedback closure. The default now uses
the original implicit conduction/reaction solver and `RecedingSlab`: retained
mass determines thickness at supplied density, the back face stays fixed, and
limited enthalpy reconstruction conservatively rezones equal-mass spatial slots.
Current geometry determines face and cell conductances. The original integrated
volatile release is transferred unchanged, while actual remap residuals enter
the combined energy audit once. The slab remains an external source to a rigid
gas vessel; its lost volume does not silently create gas volume or pressure work.

All four receding cases complete 60 seconds; the fan case also has independent
timestep and mesh refinements. Actual outputs and reproduction commands are in
`data/reference/compartment_fire/receding_pyrolyzing_compartment_example.json`.
The old fixed-geometry explicit solver and retained-cell-mass refusal remain
available through an explicit `fixed` argument; its earlier 30-second data are
preserved separately. Neither mode adds seed mass or clips release. Typed
implicit/slab refusals propagate through `VolatileGasError`, preserving direct
cancellation. Solid implicit stepping, geometry updates and gas/boundary coupling
remain first order. This homogeneous noncharring closure does not model pores,
char shrinkage, mechanics, full burnout or moving gas-boundary work.

### Finite-rate chemistry inside the opening-pressure equation

`fs-gas::multicomponent::compartment::kinetics::ZoneKinetics` can now run in
`step_with_pressure_sources`: apply the prepared operator to each transported
trial room and add its two chemical rate arrays to the callback's `additional`
account exactly once. Original adaptive homogeneous chemistry, species calorics
and common-pressure/interface work remain in their existing owners. The network
residual includes the resulting chemical pressure change; no combustion heat is
added separately. `NetworkError::Kinetics` retains typed refusals and maps
cancellation directly to the network cancellation path.

A two-room regression uses original NASA sources and a separately declared
synthetic reversible law to check finite progress, pressure-driven transfer,
element/physical-energy closure and exactly-once source accounting. It rejects
a complete chemical trial and verifies bitwise replay against a direct step
of the actual accepted shorter duration. This exercises the production owners;
the synthetic law is not methane-kinetics validation or a global entropy theorem.

### Pressure-dependent fan ventilation

`compartment::mechanical::FanModel` binds the original `FanBank`/`FanCurve` and
`LossNetwork` to their common reference density, a per-fan outlet area, total
efficiency, heat-rejection assumption and explicit density/pressure-ratio
domain. The original fan curve is static pressure. Its outlet velocity pressure
is added before matching total duct losses and actual endpoint backpressure,
including series/parallel bank scaling. Fan rise and quadratic loss scale with
actual donor density. The source convention follows ASHRAE Handbook 2016,
chapter 21, sections 5 and 10; supplied losses must include the appropriate
installation and outlet/mixing effects for that fan rating.

`ZoneNetwork::set_fan_vents` admits same-elevation point ports atomically. Each
attempt selects the actual present donor/receiver layer. Every pressure trial
recomputes the bounded operating point and transports original species and
enthalpy. Shaft input is `Q * total_pressure_rise / total_efficiency`; the
declared fraction of shaft-minus-aerodynamic loss may be rejected externally,
and the remaining power enters the receiving gas. Original enthalpy already
carries pressure-flow work, so duct dissipation is not added again. Exhaust fan
work deposited downstream into ambient does not heat an upstream room.

Only accepted operating-point diagnostics publish with room state. Nonmonotone
fan/system combinations, missing forward roots, unmodeled zero-flow power,
out-of-curve operation, domain failures and bounded-solve exhaustion refuse.
Closing a connection explicitly stops its fan and closes an ideal damper.
This is a nominal quasisteady low-pressure model with frozen inlet state during
each pressure solve, negligible endpoint kinetic storage and caller-supplied
efficiency. It does not resolve ducts, reverse flow, stall/surge or windmilling,
and retained source tolerances are not promoted to propagated uncertainty bounds.

Seven new test groups cover independent operating-point algebra, density/speed/
bank affinity, static-to-total correction, power/rejection accounting, physical
room and ambient transfers, absent-layer ports, pressure closure, cancellation
and retry. The focused mechanical filter passed eight tests including the
existing prescribed-transfer regression; full-workspace validation is separate.

### Runnable two-zone public-case comparison and pressure linearization

The `compartment_fire` example composes original NASA9 N2/O2 calorics,
Heskestad entrainment, a bidirectional stratified doorway and three original
thermal-wall profiles. Both convection coefficients are required explicit
inputs. The complete six-case Steckler comparison and dt refinement are in
`data/reference/compartment_fire/steckler_1982_two_zone_example.json`; inputs,
run command, source distinctions and all upper/lower/interface discrepancies
are documented alongside the example. Prescribed heat is not an implicit
claim of oxygen consumption or combustion-product generation.

The pressure Jacobian now shrinks its nominal forward difference with the
actual pressure residual and divides by the actual representable increment.
This prevents the fixed difference from spanning a nearby neutral-plane
transition at convergence. In the retained Test 14 scenario it reduces
whole-step halvings from 10,858 to 2 while changing the late-time mean upper
layer temperature by less than 3e-8 K. All six public runs complete with the
final solver; upper-layer-rise errors remain +4.28% to +37.79%, and interfaces
remain systematically too high. Numerical convergence is not promoted to
physical validation, and no coefficient has been fitted to these observations.

## Molecular absorption and reactive-surface radiation feedback

`radiation::molecular::MolecularPlanck` binds the original immutable gas
source and verifies the elemental/neutral identities of H2O, CO2, CO and CH4.
The corrected TNF/RADCAL Planck fits (TNF radiation archive, 7 June 2000;
Grosshandler NIST TN 1402) are evaluated on 300–2500 K without extrapolation,
coefficient fitting or silent clipping. Actual original species amounts set
partial pressures; Pa is converted to the source atm basis. Evaluation
allocates no coefficient buffers, handles absent layers and zero absorbers
exactly, and polls cancellation. A nonblank caller assumption explicitly
accounts for every unrepresented species.

These historical fits are local-temperature emission means. Using them as a
gray enclosure opacity requires a separate explicit assumption; this does not
supply spectral absorption under arbitrary incident spectra, a RADCAL solver,
a measured error bound or current spectroscopic accuracy. The original
spectral owner remains available for caller-supplied spectral data. Three
focused tests cover published coefficients, pressure scaling, amount/mass
conversion, the retained CO fit seam, domains, source identity and cancellation.

`radiation::reactive` extends the SAME existing wall/radiation feedback
iteration to original finite `Wall1d` and `RecedingSlab` owners, including
explicitly declared char density/conductivity and the existing conservative
virgin/char remap.
A slab trial advances its original implicit reaction/conduction and conservative
remap. The current remapped half-cell resistance recovers the moving-face
temperature, which closes against radiation, convection and separately
prescribed absorbed external irradiation. Only the perturbed surface response
is recomputed for a Jacobian column; original prepared workspaces are reused.

An accepted candidate retains the same signed query powers in every condensed
owner and opposite gas account. Endpoint face-temperature, optical-power,
total-boundary-power and combined physical-energy residuals are checked
separately. Original volatile mass and enthalpy remain outgoing in the slab
report and must enter the gas through `VolatileMixture` exactly once;
`report.add_gas_heat` contains only radiation and convection. Install all
returned owners only with their accepted network interval, and warm-start from
`candidate.input_powers`, which contains TOTAL slab front heat.

Four focused reactive cases cover original inert-wall agreement, reactive
mass/enthalpy/flux closure, actual NASA volatile transfer with network retry and
deterministic replay, spectral coupling, work limits and cancellation. The
native airflow library passes all 110 tests. Slab backs remain insulated;
gas state/opacity and optical geometry are frozen over each interval. Recession
does not infer gas-volume change or pressure work. This is first-order split
coupling, not a guaranteed second-order or fully implicit gas integrator,
a coupled entropy theorem, flame structure, soot production or experimental
compartment-fire validation.

## Coupled full-methane chemistry and moving retained char

The existing `coupled_compartment_fire` example now selects `hydrogen`,
`methane` or `charring-methane` through one fuel/source adapter. The latter
two use the complete pinned 53-species, 325-channel GRI-Mech 3.0 mechanism,
including its original irreversible channels. Canonical species indices,
molar masses and element inventories come from the admitted source, and all
species stay in the original gas/network owners. Selected CO, CO2 and NO are
reported without removing other intermediates or inferring soot yields.

Methane starts in explicitly hot 1500 K oxygen/nitrogen gas with no radical
seed. The gas admits 300--3000 K; the existing molecular Planck closure still
enforces its narrower 300--2500 K range. CH4, CO, CO2 and H2O contribute their
actual partial pressures; the unrepresented species remain explicitly
transparent. Optical state is frozen per outer interval as before.

The manufactured charring fuel is defined by mass as 0.8 CH4 plus 0.2 condensed
carbon. Its independent char enthalpy, 300 kg/m3 char density, 0.08 W/(m K)
char conductivity and parallel volume-fraction conductivity mixing feed the
same conservative `RecedingSlab::with_char` owner and implicit radiating-face
solve. The initial virgin enthalpy is the yield-weighted original volatile/char
enthalpy at 300 K minus the authored 1 MJ/kg decomposition energy. Actual
time-integrated volatile mass and enthalpy enter `VolatileMixture` exactly
once; source chemical formation energy is never added as separate combustion
heat. Independent audits include gas, remaining virgin, retained char, finite
wall storage and external vent streams.

The new modes share the original pressure-source callbacks, prepared
`ZoneKinetics`, nonlinear surface coupling and accepted-payload publication.
No thermal, chemistry or transport clock advances twice. The summary includes
accepted chemical work and checked direct-Newton fallback counts. Counts from
discarded pressure-root trials are outside that accepted-payload summary.

The existing hydrogen SDIRK2 consumer and both frozen charring controls were
executed through 0.1 s with 0.001 s maximum step and 12 solid cells. Hydrogen
reproduces its prior 1331.048718306 K, 137365.531842355 Pa endpoint. Frozen
closed and vented controls each complete 100 outer intervals and conservative
remaps, with no network halving or gas reaction products. Their maximum combined
residuals are below 6.43e-17 kg mass, 3.26e-7 J energy and 1.41e-14 mol elements.
Both closed reacting methane variants also complete 0.1 s in 100 outer
intervals with no network halving. Their final temperatures are 1683.198402902 K
without char and 1646.380168650 K with char; the latter retains
3.691197223245e-6 kg of condensed carbon. The largest combined residual across
these closed reacting cases is 7.16e-17 kg mass, 3.26e-7 J energy and
1.20e-14 mol elements. A separately completed 20 ms reacting vent case uses
20 outer intervals with no halving, reaching 1522.866245311 K at
102199.510551066 Pa; its maximum pressure residual is below 8.51e-6 Pa.

Commands, source inputs and actual endpoint tables are retained in
`examples/coupled_compartment_fire/README.md`.

This is a manufactured source-composition exercise. It is not a measured wood,
PMMA, hydride, char-oxidation or pyrolysate mechanism. Recession changes only
the solid geometry and supplies no inferred vessel-volume displacement or
pressure work. The existing first-order overall source/pressure/thermal
composition, approximate gray radiation, finite chemistry domains and absence
of an experimental fire/soot/emission validation remain explicit.

### Original porous permeability in liquid operating points

The native temperature-feedback pressure path now accepts a sealed resistance-owner
abstraction while preserving the existing geometry-law default. Original continuum
permeability supplies the unchanged tensor/series sum of thickness over effective
permeability divided by flow area. An explicitly calibrated inertial skeleton
also supplies `sum(beta L)/area²`; no beta is inferred from permeability or pore
radius. Both use the original `dp = a mu Q + b rho Q |Q|` hydraulic owner with
the same current fluid snapshot. That owner evaluates current
liquid density, viscosity, Reynolds admission, pressure loss and dissipated power.
Mixed ordinary/porous sections, the fully dissipative enthalpy march, coupled-solid
pressure closure, parallel branches, pump-coupled networks, recirculating coolant
and stored solid/contact/radiation transactions share that owner.

A porous binding retains the original layers, frame, source claims and support and
requires an explicit absolute downstream path pressure. Every local hydraulic
sample admits its actual bulk temperature and superficial speed. Every complete
trial (including bracket endpoints and rejected solid/flow candidates) reconstructs
the actual regional absolute pressure faces by a reverse stream walk and admits
those faces, fluid endpoint temperatures/speeds and bulk state against the original
support. Porous sections and parallel manifolds must agree on the outlet datum.
Accepted regional reports retain the original barrier and actual pressure faces.
Changing permeability or liquid transport data therefore propagates through flow,
heat transport, solid feedback and enthalpy-conserving outlet mixing.

Scope remains the native positive-flow, low-compressibility regional liquid model:
transport pressure stays at its declared reference, intrinsic permeability is
constant within its source support, and the skeleton is rigid and fully liquid
filled. Gas slip/relative saturation, a sealed positive-flow path, inconsistent
pressure datums and unsupported trial states refuse without coefficient floors.
No liquid wetting/multiphase mobility, evolving pore geometry, skeleton heat store,
pressure EOS or inertial porous law is inferred. The opt-in hydraulic-heating path
retains its explicit complete conversion to fluid heat; the other paths retain
their original omitted-power allowance.

Regression sources exercise data-only permeability/viscosity changes, independently
derived pressure/heat operating points, hydraulic heat in the same enthalpy state,
original source receipts, discarded-trial support refusal, sealed topology, solid
feedback and parallel redistribution/mixing. The calibrated-inertia regressions
independently reconstruct liquid operating points and regional exponential
enthalpy balances, including temperature-dependent viscosity and the inertial
heat contribution exactly once. Gas graph regressions exercise the corresponding
signed flow, fan intersection, coefficient binding and actual source limits.


## Complete coolant assemblies on native solid motion

The existing `SourceConvectionHeatNetwork` and its finite `with_links` view
also serve the shared `observe_heat_assembly` timeline and all four native
static/dynamic thermoplastic and reciprocal thermoelastic assembly entries.
Source-bound coolant and original solid material temperatures remain on the
same native candidate. Finite coolant storage/transport can therefore change
actual deformation and retained plastic history while preserving original
material identities, heat destinations and source receipts.

The complete source sample owns the explicit additional-rate bound: original
radiation plus prescribed reservoir h*A, with native finite-link conductances
excluded because they already enter the thermal owner. The specialized
convection observer and generic assembly observer reuse this same calculation.
The independent combined-source/event regression checks identical substeps,
continued history, finite fluid endpoints and external heat; an undeclared
explicit bound refuses without changing the original history.

The finite-coolant J2 regression changes only original coolant Cp data and
compares actual loaded equilibrium with an independent two-shear return
mapping, then continues cooling through unloading with retained plastic
history. Physical heat, finite internal exchange, original mechanical storage
and numerical loss remain separate. Fixed geometry and the existing thermal
and mechanical refinement requirements continue to apply.

### Prescribed liquid saturation in porous pressure/heat paths

`PorousLiquidResistance::with_liquid_states` updates each original porous layer
with its explicit wetting-liquid calibration through the shared continuum
geometry. The native pressure/enthalpy solver uses sum(L/(k_effective*kr_l))
and the same current liquid density and viscosity already used by its thermal
operating point. Flow area, Reynolds length/range and the common absolute outlet
datum are retained. Repeated updates replace the old saturation states and
recompile from intrinsic k, never compound earlier mobility reductions.

Every accepted porous regional report retains the original barrier plus the
exact ordered `liquid_states` used by that trial. Actual bulk and endpoint
temperatures, both absolute pressure faces and superficial speeds must satisfy
the intersected original-skeleton and liquid-calibration limits. Unsupported
bracket trials refuse even if a proposed root would be admissible. Exact dry
seals refuse positive-flow series paths; no permeability floor opens them.

Capillary pressure remains Pc=pg-pl from the retention law. Changing its
calibration at prescribed saturation cannot create another hydraulic pressure
loss or thermal source. The existing declared dissipation convention remains
responsible for pressure-work/heat accounting. This is an instantaneous
operating point on prescribed saturation: liquid storage, gas pressure,
Richards evolution, wetting fronts and hysteresis are not advanced.

### Spatial material response within native graph branches

`natural_circulation::distributed_response` compiles independently sampled
hydraulic sections on each original branch. Original bulk viscosity supplies
`A = sum(a_i mu_i)`; the conserved reference density supplies
`B = sum(b_i rho0)`. The native graph then uses the same signed law
`dp = A Q + B Q |Q|`. A cold, viscous section therefore affects the actual
network flow instead of disappearing into one branch-average property.

The adapter reuses the existing material/correlation compiler once, with
bounded linear section work and compensated positive resistance sums. The
temporary batch graph is never solved. Original nodes, endpoints, branch
identities and declared exchange order are preserved. Every section retains
its material owner, original law, source classification, loss evaluation and
heat evidence. Native summed coefficients carry an explicitly derived
engineering allowance; the original evidence is not promoted.

Authored section factors and wetted areas must partition the actual geometry.
Subdivision never restarts a developing-flow inlet. Material samples, spatial
refinement and integrated gravity remain the coupled owner's responsibility.
The original source/Reynolds/correlation restrictions apply at every query;
unsupported states and over-budget section counts refuse. Regression sources
cover unequal-viscosity parallel flow, reversal, geometry subdivision, actual
exponential wall heat and convergence toward independently integrated drag.

### Retained momentum in branched incompressible graphs

`graph::inertial::InertialGraph` retains one branch impulse `L Q` with
`L = rho0 length / area` and kinetic energy `L Q^2 / 2`. Native
`DescriptorPortHamiltonian` incidence stamps enforce free-node continuity;
pressure multipliers have zero storage. Multiple independent circulation modes,
pressure ports and additional signed branch heads share the existing native
descriptor step. No pressure penalty or separate momentum integrator is added.

Immutable proposals may resolve new original linear/quadratic coefficients and
Boussinesq density samples while preserving branch identities, reference
density, inertia, pressure-node identities and gravity geometry. Native reduced
pressures retain hydrostatic rest at uniform density. Dedicated reports expose
interval pressures, actual midpoint flows and final retained flows separately.
Average rates are the transport interface; final rates retain momentum history.

Acceptance independently gates the original nonlinear branch impulse residual
and mechanical energy balance. Continuity, cancellation, native Darcy support,
bounded descriptor dimension and per-branch drag-step limits remain explicit.
Non-cancelling omitted-work bounds use the largest old/new/average speed, so a
zero average during reversal cannot erase unresolved dissipation. Broader source,
phase, Reynolds and correlation admission remains with the original material
callback. Pressure work, kinetic storage, dissipative loss and anomalous-gravity
work remain separate accounts.

The current native dense owner is limited to 256 total branch/multiplier states.
It is an incompressible, fixed-reference-density model, with no inferred acoustic,
compressible, turbulent or multiphase storage. Native regression sources cover
branched startup/decay/reversal, nonlinear loss closure, source gas buoyancy and
viscosity, explicit pressure-port work and transactional refusal.

### Pipe support forces and moments from native momentum

`InertialGraphTrial::pipe_load` also converts that same trajectory to a physical
support force, torque and integrated impulses for each original fixed pipe.
`PipeLoadGeometry` compiles the actual constant-section centreline once; its
length and area must independently match the original inertia. Bound gravity
also requires the original endpoints and every distributed source section.
The reaction includes individual endpoint pressure tractions, accelerating
retained fluid, advected momentum and original local material weight. The
latest constitutive recheck supplies gravity samples. A declared uniform
exterior pressure shares the nodal datum, so shifting a gauge cannot create
force. Full line and weight moments retain bends and arbitrary moment origins.

Momentum flux integrates the square of the native affine flow, including both
directions of a reversal with zero average transport. Friction is already in
the solved pressure trajectory and is not added again as a second wall force.
The returned wrench is directly usable by a mechanical load/impulse owner;
regressions include an original articulated-body impulse and independent work
balance. This is a slender plug-flow load on a fixed control volume, with no
inferred shaft torque, swirl, reducer/junction closure or moving-pipe feedback.
Publish candidate loads only after that same hydraulic trial is accepted.

### Multiport fitting and complete hydraulic assembly loads

`InertialGraphTrial::junction_load` now resolves each original zero-storage
junction's force, torque and impulses from every incident pipe. A starting
pipe's outward junction normal is its first tangent and its outward flow is Q;
an ending pipe uses minus its last tangent and -Q. The endpoint positions must
coincide, the original area and length remain admitted, and both endpoint flow
balances must close. No additional junction mass, gravity load or species/
enthalpy transport is created. The common STATIC pressure remains the original
graph approximation; unequal areas do not imply ideal reducer recovery.

Each junction retains signed outward kinetic-energy transport
`integral sum(rho*q_out^3/(2*A^2)) dt`, using the exact cubic integral of the
same affine native flow trajectory. The quadratic momentum integral remains
nonzero through a zero-average reversal. Signed energy transport is an
unresolved physical term, not a passive-loss claim or an automatic heat credit.
Positive and negative terms require their own pressure/thermal closure; equal-
speed turning can close this term naturally. Per-junction values remain visible
even when their assembly sum cancels.

`assembly_loads` retains original branch-order pipe loads and node-order fitting
loads. Its total omits internal endpoint tractions/fluxes analytically before
summation, retaining all stored-fluid acceleration and original local weight.
It can therefore recover small net support loads without subtracting large
internal pressure forces. Apply either the separate component loads or their
total, not both. Free nodes are internal; pressure nodes remain external stream
cuts unless explicitly declared enclosed pressure datums. Such a declaration
means no external stream and independently requires endpoint continuity.

`accept_with_assembly_loads` publishes the hydraulic state and load package only
after the original impulse, mechanical-energy and cancellation gates pass.
The same fixed-geometry and one-way mechanical scope as pipe loads applies.
Native regression kernels exercise original descriptor startup and an actual
articulated support impulse, equal-speed turning, unequal-area signed energy,
reversal/subdivision, closed-loop fluid spin recoil and exterior-load cancellation.

### Joint branched momentum and finite enthalpy convection

With `transient-transport`,
`thermal::transport::transient::inertial::advance_inertial_graph` couples the
retained branched momentum owner to the original finite-cell backward-Euler
thermal owner. The native Newton--Krylov residual closes cell energy increments
and candidate-versus-native branch-average rates together. A pure constitutive
callback resolves actual candidate material states; the actual-flow callback
re-admits retained sources at native averages and old/new speed envelopes,
then supplies wall exchange and explicit external mass/enthalpy streams.

`GraphEnthalpyTransport::from_volume_flows` binds those actual average rates
directly to graph topology, preserving the original reference inventory,
caloric identity and all-cell mass closure. No steady pressure solution is
fabricated for accelerating fluid. Internal upwind transfers remain internal;
external advective heat and native nonadvective heat enter exactly once.
Unmapped solid cells and their original conduction/boundaries remain active.

A fresh final proposal independently gates flow closure, original impulse,
mechanical energy and per-cell/total thermal energy. Both original inputs are
immutable; updated momentum and thermal states return at one accepted time
only after every gate and final cancellation checkpoint passes. Rest preserves
finite heat storage. A reversed endpoint velocity does not reverse net heat
transport when the native interval-average velocity is still positive.

Midpoint momentum with endpoint-held backward-Euler heat is first order as a
coupled method. Time and cell refinement, source support, original thermal
event splitting, explicit low-Eckert caps and bounded native dense momentum
remain required. No reversal-time event locator, compressible energy equation
or inferred spatial mesh is added. Five regression sources exercise independent
branched startup and thermal roots, directional heat transport through reversal,
temperature-dependent drag feedback, rest/refusal behavior and explicit
reservoir mass/enthalpy binding.

### Actual-flow constitutive closure for inertial convection

Final inertial thermal acceptance refreshes original graph and source
coefficients at the actual native interval-average rates. A small
candidate-versus-actual flow difference alone cannot establish closure of a
sharp supported source law. `InertialGraphTrial::recheck` audits the exact
retained impulse trajectory against the refreshed loss, body-head and pressure
data without advancing momentum again. The original impulse, work, support and
speed-envelope gates are shared with ordinary proposal evaluation.

The refreshed source payload also drives a fresh native thermal evaluation at
the same enthalpy field. This evaluation is charged to the existing work budget.
The native descriptor receipt still describes its original held-coefficient
solve; the independent original-law ledger describes the refreshed accepted
physics. Regression sources distinguish these accounts and reject a sharp
continuous drag law whose rate mismatch is tiny but whose impulse and energy
defects are significant. A failed refresh publishes neither physical clock.

### Sparse frozen-flow implicit enthalpy transport

`GraphEnthalpyTransport::advance_advection_sparse_implicit` retains the original
native backward-Euler enthalpy residual while assembling its exact upwind
energy-increment operator with native COO/CSR storage. Native flexible GMRES
with a diagonal preconditioner replaces the small dense solve for larger mapped
inventories. Matrix entries include the native component-reference derivative,
so large formation enthalpies do not become artificial heat sources. Reverse
flows and open mass/enthalpy ports preserve the original donor orientation.

Explicit limits bound cells, staged sparse entries, Krylov workspace, cycles and
representability corrections. Native per-cell and total joule residuals remain
the final acceptance gates; a linear-solver receipt alone cannot advance the
thermal state. This path supports held mass rates and fixed cell sources, with
all active cells mapped. Active conduction, Robin exchange and reservoir laws
use the existing general endpoint interface instead. Implicit stability does
not remove timestep or spatial-refinement requirements. Regression sources
exercise 96-cell open and closed inventories, nonlinear caloric inversion,
reversal, formation-datum shifts and transactional work-limit refusal.

### Distributed pipe residence on coupled graph timelines

`transient::distributed::DistributedGraphInventory` binds ordered original
thermal control volumes to each hydraulic branch. Only native mass-transfer
connections expand; junction pressure and branch momentum dimensions stay fixed.
Every pipe cell retains its own enthalpy through changing and reversed flow.
Explicit cell geometry supplies volume, area and the equivalent branch inertia
`rho0 sum(length_i / area_i)`, including nonuniform cross-sections. Inventoried
branches must agree with the original momentum owner; empty paths supply no
inferred volume or inertia. External streams enter only declared junctions.

The distributed inertial entrypoint uses this same inventory
at every thermal candidate and final source refresh. Original node order, branch
identity, reference density, full cell geometry, mass and available native caloric
identity remain bound. Returned transport bindings also reject replaced geometry
or caloric identity. The native thermal owner still owns all cell energies,
wall links, events and acceptance gates. No second time integrator is introduced.

This is a finite-volume mixed-cell axial model. The caller authors its spatial
resolution; nominal V/|Q| residence is a held-flow diagnostic, not a prediction
of breakthrough under changing flow. Time and spatial refinement remain necessary.
The distributed inventory admits at most 4096 fluid cells; original nonlinear,
momentum and thermal work limits apply independently. Regression sources exercise
pipe residence and spatial refinement, variable-area inertia, retained heat during
reversal and rest, coupled timeline acceptance and transactional binding refusal.

### Distributed storage with quasi-steady hydraulic feedback

`advance_quasi_steady_graph_distributed` uses the same distributed pipe inventory
on each original backward-Euler thermal candidate. The existing quasi-steady
entrypoint and problem structure remain compatible. No hydraulic nodes or pipe
cells are inferred. Actual solved signed branch rates are expanded through the
authored spatial order, including reverse flushing of already heated fluid.

Node mapping, reference density, original topology, native storage identity and
expanded connection budgets are admitted before source callbacks. Each actual
solved gravity chart is checked against the inventory as well as the original
graph template. Native thermal energy, source support and work budgets retain
the existing acceptance rules. Hydraulic momentum is still quasi-steady; the
inertial entrypoint owns retained branch impulses. Explicit distributed gravity
retains absolute local anomalous heads before branch cancellation; the native
work ledger uses that sum. Original source callbacks independently admit section
drag, gas Mach and local mechanical work at the actual solved rates.

Regression sources couple temperature-dependent source viscosity to finite pipe
residence and reverse flushing, checking an independent ordered cell recurrence
and scalar flow root. Geometry, caloric identity, node/density/topology, actual
gravity and connection-budget refusals preserve the original thermal timeline.

### Original liquid and gas sources on distributed graph storage

`transient::material::{LiquidGraphSources, GasGraphSources}` bind original source
laws to the exact native caloric object registered for every fluid cell. Each
junction and pipe temperature is resolved through the existing material owner.
Liquid pressure and source coordinates remain intact; gas NASA enthalpy,
formation datum, thermodynamic pressure and composition remain intact. Reference
rho0 owns conserved mass, continuity and inertia. Local viscosity owns drag;
local density departure owns Boussinesq gravity. This remains a bounded low-Mach,
low-Eckert approximation with explicit original source support.

Each hydraulic branch has one original adiabatic resistance law per pipe cell.
The source adapter retains every actual pipe-cell start/end and material sample.
Bent branches use the sum of `(rho_i-rho0) g.dx_i`; equal-height endpoints and
cancelled diagnostic mean density cannot erase a hot/cold circulation head.
Length-weighted branch densities are diagnostics only. Construct retained
momentum from the initially resolved source graph,
using the distributed inventory's equivalent geometry, so its reference density,
gravity chart and admission cap are already fixed before the first coupled step.
The source payload belongs in the inertial driver's resolve result; its exchange
callback calls `admit_actual` or `evaluate_links` before returning coefficients.

The same source owners expose `admit_quasi_steady` and
`evaluate_quasi_steady_links` for an actual native `GraphSolution` and its
positive thermal interval duration. A shared borrowed flow view retains the
native solution instead of manufacturing an inertial report. Quasi-steady rates
are held over that thermal interval; local power times that duration supplies
the omitted-energy bound. Original topology, gravity frame, exact pipe path,
reference mass, source identity, Reynolds support and gas Mach remain admitted.

Actual-flow admission rechecks original section laws at native old, new and
average rates, including zero-crossing support. Maximum speed independently
bounds gas Mach and omitted mechanical work. Absolute local anomalous-gravity
work is summed before sections can cancel, then combined with hydraulic
dissipation under separate power and energy caps. Prepared source payloads are
tied to the exact binding and candidate temperature field. The coupled driver's
final actual-flow refresh retains these original source checks.

`material::convection` reuses original reciprocal finite-cell heat links. Forced
speed is the actual reference-volume rate divided by authored native pipe area.
Bulk/film rho, mu, k and Cp come from the same original source and correlation.
Support is re-admitted at both endpoint speeds and any signed crossing; an
admissible average or maximum alone cannot establish a lower Reynolds limit.
Samples are associated by native link identity regardless of caller ordering.
Pure natural cards require a quiescent bulk-flow envelope. The explicit
`GraphConvectionLaw::AidingVerticalPlate` uses the existing Manning-Qureshi
unconfined laminar relation, actual signed cell-tangent velocity, original film
properties and sourced thermal expansion. `AidingVerticalPlateFrame` admits the
authored upward plate axis against actual gravity and relative flow; inclined
plates, crossflow and opposing buoyancy cannot borrow the aiding relation.
Pure limits and resting thermal equilibrium reuse the existing finite-cell
owner. Old/new rates and zero crossings must satisfy that same original law.
Apply returned conductances to the native link map once, without duplicating
their heat as cell sources.

Regression sources cover spatial viscosity, original EOS/calorics, source and
payload identity, non-cancelling local work, endpoint Reynolds/Mach support,
reciprocal wall exchange and native coupled operation. No new experimental
coverage or material/correlation claim is inferred by the adapter.

### Distributed gravity on original hydraulic and momentum owners

`LossGraph::with_distributed_boussinesq_gravity` retains a nonempty ordered
`GraphGravitySection` path for each original branch. Sections must be finite,
continuous, positively sized and terminate at the original graph nodes. Every
local density, including horizontal sections, obeys the explicit Boussinesq cap.
An explicit total section budget is retained with the path. Source refresh may
change material samples but cannot change a retained momentum owner's path,
partition or section budget. Distributed thermal inventory additionally binds
one gravity section to each stored pipe cell's exact geometry.

Binding caches both signed anomalous line integrals and sums of absolute local
anomalous heads. Existing steady and inertial solvers consume those caches;
there is no added hydraulic mesh, invented effective density, external drive or
replacement momentum solver. The existing hydrostatic reference frame preserves
uniform-density rest. Absolute local heads feed both native quasi-steady and
inertial omitted-work gates, even when their signed branch sums vanish.

### Bounded native time resolution through flow reversal

`transient::inertial::adaptive::advance_inertial_graph_resolved` compares one
original coupled step with two native half steps, then retains only the two
independently admitted fine steps. It supports original junction storage and
explicit distributed pipe inventory. Every active cell energy and every branch
flow participates in the local estimate. Relative thermal tolerances scale
energy changes, preserving independence from the formation-energy datum.

Accepted substeps also obey explicit limits on sign-changing interval duration
and the cancelling counterflow volume of their affine momentum trajectory.
Thus a nearly zero average rate cannot silently erase both directions of
transport. Native/source callback budgets count discarded comparisons as well
as retained steps. Source, correlation and physical refusals remain terminal;
the entire requested window returns synchronized owners or leaves both input
owners unchanged, including after a locally accepted prefix.

These controls estimate local discretization error and resolve sampled reversal
trajectories; they do not certify global error or locate events exactly. The
caller still authors spatial resolution, maximum interval and scheduled-event
splits. All thermal conservation, source support and momentum work gates remain
on the original native reports.

### Source-backed axial conduction in distributed fluid storage

Liquid and gas graph source owners expose `evaluate_conduction_links` for
explicit native zero-baseline links between consecutive cells of the same
original branch. Each cell's own bulk source conductivity and actual half length
and cross-sectional area determine the series resistance
`L_a/(2 k_a A_a) + L_b/(2 k_b A_b)`. The adapter reuses the original thermal
series-conductance kernel and reciprocal link integration. Conduction remains
active at rest and independent of advective direction; diagnostic signed heat
must not also be supplied as an external cell source.

The declared face must cover the smaller adjacent whole pipe section. Partial
apertures, junctions, cross-branch links and duplicate whole-face owners refuse.
Stepped sections use a one-dimensional held-cell model with no inferred
spreading, constriction, contact resistance or turbulent dispersion. Exact
candidate-temperature and original source-binding admission prevent stale
coefficients. Immutable source sections are shared with the accepted payload
without copying source receipts per face. Explicit link budgets bound selected
links and the linear scan of active native links; cancellation and conflicting
coefficient ownership leave the destination map unchanged.

Regression sources cover variable conductivity and unequal sections, independent
Fourier-profile mesh refinement, stagnant native equilibration, signed advection
with conduction, original nonlinear gas calorics, and geometry/ownership refusal.
Time and spatial resolution remain caller responsibilities.

### Limited heat-front transport on the shared graph endpoint

`DistributedGraphInventory::with_limited_reconstruction` prepares native
minmod enthalpy traces only for explicitly ordered interior pipe cells. Actual
adjacent lengths determine the slopes, including nonuniform and bent paths.
Junctions and branch endpoints retain original donor-cell traces. Signed solved
flow selects the actual donor face, while one conservative native transfer owns
both cell debits. The immutable geometry plan is shared across bindings without
copying it on every coupled evaluation.

Existing distributed quasi-steady, inertial and resolved-time drivers receive
the same reconstruction through `GraphEnthalpyTransport::endpoint_transport`.
`GraphEnthalpyTransport::advance_implicit_with_heat` additionally solves prescribed
signed flow with endpoint-resolved wall/source coefficients through the existing
Newton/Krylov owner and original backward-Euler energy residual. Its bounded
endpoint-query count includes numerical derivatives and final physical
acceptance. Heat callbacks supply nonadvective coefficients; native internal
transfers and external stream enthalpy are included once.

The same authored plan now reaches the native explicit timeline through
`GraphEnthalpyTransport::wrap_observer` and `advance_to`. The shared
`ReconstructedAdvected` observer applies the original signed face correction
while retaining source/conduction work, scheduled events, native calorics and
the combined convection/conduction step restriction. No extra heat integrator
or per-step reconstruction-plan clone is introduced. A missing or empty plan
uses the original upwind arithmetic. Nonuniform affine transport with original
external heating, both flow signs and rest has a direct regression source.

A nonempty reconstruction still cannot enter the linear implicit shortcut:
its upwind Jacobian does not describe the nonlinear limited flux. Storage,
original caloric identities, source support, scheduled-event splits and
independent local/global joule acceptance remain authoritative. The original
state and payload remain unpublished on failure, source refusal, cancellation
or exhausted work.

Regression sources compare smooth-front transport against an independent
continuous plug-flow profile, connect prescribed and coupled native owners, and
check reversal, rest, datum shifts and transactional refusal. The limiter bounds
individual face traces; it does not establish a global maximum principle for
arbitrary coupled sources or unstructured junctions. Both the native explicit
stepper and backward Euler remain first order in time. Spatial/time refinement, including reversal resolution,
remains necessary; no measured throughput or unconditional nonlinear convergence
is claimed.

### Prescribed pump and inlet histories on graph thermal transport

`GraphEnthalpyTransport::with_schedule` binds the original thermal
`AdvectionSchedule` to the admitted signed graph and fixed fluid inventory.
The shared history scales every original internal and external mass rate by
one nonnegative multiplier and samples inlet enthalpy in its original caloric
law. Baseline rates remain immutable, so stopping the pump does not erase the
flow used on restart. Source, material, inlet identity and storage admission
remain with the existing native owners.

`wrap_observer` and `advance_to` carry that history into the original explicit
timeline, with or without limited axial reconstruction. Each accepted substep
uses its start-time sample, the actual sampled CFL restriction, and the earliest
native source, schedule or authored-event transition. Existing dense, sparse
and nonlinear implicit graph entries instead freeze one endpoint sample before
their native backward-Euler solve. Newton/Krylov trials reuse that one flow;
they do not repeatedly sample or integrate the schedule. An implicit interval
containing an interior schedule knot refuses before advancing either owner;
callers split at the knot. An interval ending exactly at the knot is admitted.

Native reciprocal transfers, external stream enthalpy and nonadvective heat
retain their original single owners, conservation checks and rollback.
Unscheduled bindings preserve their existing execution paths. Regression
sources cover a stopped/ramped two-tank analytic response and refinement,
agreement of all three implicit entries with a common endpoint balance,
datum shifts, invalid inlet identities, interior-knot refusal, and both flow
signs in a reconstructed distributed pipe.

The schedule prescribes flow; it does not solve pump momentum, infer mechanical
work, reverse branches independently, change stored mass or evolve composition.
Explicit start sampling and backward-Euler endpoint sampling retain their
original first-order time accuracy. Callers still choose and refine intervals
that resolve pump startup and inlet variation. No measured throughput or global
temporal-error bound is implied.

## Actual liquid graph dissipation as local fluid heat

`LiquidGraphSources::dissipative_heating` binds an opt-in heat owner to the
original distributed pipe inventory and liquid caloric source. The original
constant-density dissipation admission qualifies the whole caloric source
span. Actual node and section source densities must still equal the retained
reference density at evaluation. Variable-density liquids and gas sources do
not inherit this pressure-work closure.

`LiquidGraphHeating::evaluate_inertial` and `evaluate_quasi_steady` first use
the existing source/flow admission, retaining source identities, exact pipe-cell
paths, original resistance support at old/new/average and crossing flow, and
all local source and native mechanical-work limits. Each original section's
actual dissipated power is deposited in its own fluid cell. The source graph's
conservative maximum-speed envelope remains an applicability bound; it is not
the heat supplied to the fluid.

The resolved receipt retains original source sections by shared ownership,
per-section watts and joules, the native source admission, mechanical debit,
thermal credit and transfer residual. Both the native sum of branch pressure
drop times actual average flow and the local heat sum must agree with the
inertial report's hydraulic-loss debit within a roundoff bound. Quasi-steady
evaluation uses the corresponding original solved-branch power over the
declared interval. A changed debit, stale source endpoint, nonfinite or
unrepresentable transfer refuses the entire receipt.

`ResolvedLiquidGraphHeating::apply_to` adds this field to the native additional
cell-source map only after checking every ownership conflict and cancellation.
Zero-power samples remain diagnostic but claim no additional-source slot.
Original authored cell heat is retained by the native thermal owner. Callers
can compose the receipt in the existing inertial or quasi-steady thermal
exchange callback, so original enthalpy storage, finite pipe residence,
advection and reciprocal wall heat share the same solve and accepted clock.
No second integrator or synthetic resistance is introduced.

Each modeled hydraulic section is assumed to deposit its entire irreversible
loss into its represented liquid cell. A porous skeleton or separate body
receiving a share requires an explicit partition model. Kinetic storage,
pressure-port supply, anomalous-gravity work, pump inefficiency and numerical
dissipation are not thermal sources in this adapter. Existing low-Eckert and
source/native work ceilings remain conservative admission limits.

Regression sources cover independent local losses with viscosity variation and
flow reversal; an unforced closed loop's coupled kinetic-to-thermal transfer;
and source-density, stale-epoch, duplicate-owner, native-debit and work-limit
refusal. Compilation and native numerical execution were unavailable in the
delivery environment.

## Prescribed-wall convection on original graph fluid storage

`material::convection::reservoir` supplies source-backed conductance for an
already-present native reservoir on a bound liquid or gas cell. Both source
owners expose `evaluate_reservoirs` for an actual inertial report and
`evaluate_quasi_steady_reservoirs` for a solved graph and its thermal interval.
The reservoir retains its declared covered area and prescribed temperature.
Its authored conductance must be zero; missing reservoirs, duplicate bindings
and competing conductance owners refuse without changing native state.

Bulk sampling uses the actual finite fluid temperature; mean-film sampling
uses that temperature and the actual wall temperature. The full wall-to-fluid
temperature difference drives heat and buoyancy. Original liquid source laws
or fixed-pressure/composition gas sources retain their transport properties,
thermal expansion, support and registered caloric identity. Source samples,
native geometry, candidate temperatures, density limits and original mechanical
work gates remain admitted. Gas evaluation retains actual endpoint Mach limits.

Finite-wall and prescribed-wall convection share the original flow-law
admission. Forced cards use actual signed reference-volume flow and pipe area;
old/new/average rates and any zero crossing retain original regime support.
Pure natural cards require a quiescent original flow envelope. Aiding mixed
plates retain original local tangent, gravity orientation, film properties and
density-departure restrictions. No new correlation, empirical blend, velocity
floor or substitute material tuple is introduced.

`ResolvedReservoirConvection::apply_to` supplies total W/K values to the native
`BoundaryConductances::reservoir_w_k` map. It checks all ownership conflicts and
cancellation before mutation. The existing thermal boundary owner accounts for
heat, external energy and entropy once; diagnostic sample watts must not also
be added as cell sources. Existing ambient patches, authored sources, finite
contacts, native enthalpy storage and scheduled boundary events keep their
owners. Re-evaluate conductance at each coupled thermal trial and apply it to
that trial's coefficient map. No second thermal integrator is introduced.

Explicit reservoir counts bound work and storage; samples are ordered by native
cell identity. Nonfinite or unrepresentable transfer refuses. Four regression
sources cover independent forced conductance and full wall/bulk temperature
drop under heating/cooling and signed reversal; a native coupled natural-
convection endpoint against an independent caloric balance; boundary ownership,
stale source and zero-crossing refusal; and original gas film EOS, NASA calorics
and endpoint Mach admission. Source/API review was completed. Rust parsing,
compilation and native numerical execution were unavailable in the delivery
environment, so these regression sources are not reported as executed.

This is exchange with a prescribed external wall temperature. Finite wall/body
storage remains on original reciprocal links. Original material/correlation
applicability, low-Mach/Boussinesq assumptions and time/spatial resolution remain
caller responsibilities; no new experimental calibration or throughput claim
is made.

## Developing convection on native graph sections

Liquid and gas graph source owners expose `evaluate_developing_sections` and
`evaluate_quasi_steady_developing_sections`. The shared
`material::convection::developing` implementation uses the original Hausen mean
for `DevelopingDuct::Circular`, or the Shah-London simultaneously developing
rectangular CWT table for `DevelopingDuct::Rectangular`, through `fs-convection`'s
section-integrated evaluators. Explicit forward/reverse
entrance distances follow actual signed flow; spatial subdivision does not
implicitly restart the heated entrance at every pipe cell.

Each binding names an original pipe cell and an existing zero-conductance
finite wall link or prescribed wall reservoir. The authored diameter or
width/height must match the native flow area, and the wall patch must cover
that cell's full wetted perimeter and length. Rectangular geometry uses
`Dh=4A/P` in both Re and `h=Nu k/Dh`, and the physical minor/major side ratio
in the table admission. One binding owns that whole wall.
Original local bulk/film sources, caloric identity, density, Reynolds/Mach and
mechanical-work limits remain authoritative. Old/new/average flow and signed
zero crossings retain original correlation support and entrance orientation;
the forced card still refuses rest and trajectories crossing zero Reynolds.

`ResolvedDevelopingConvection::apply_to` supplies native link and reservoir
conductances after checking both destinations for conflicts and cancellation.
The existing thermal owner performs reciprocal wall transfer, external
reservoir heat, advection and enthalpy storage once on its original clock.
Returned wall-to-fluid powers are diagnostic and must not be supplied again
as cell sources. Explicit section counts bound source work; native cell order
fixes result ordering. No new momentum or thermal integrator is introduced.

The circular relation retains its constant-wall-temperature, hydrodynamically
developed inlet. The rectangular relation instead retains simultaneously
developing flow, equal constant temperature on all four walls, exactly
aspect ratio 0.5 and the actual source Pr=0.72. Unsupported properties are
refused, not rounded or replaced with a circular coefficient. Original table
and declared low-Graetz engineering-bridge evidence stays on the retained
`DevelopingSectionEvaluation`. The caller retains the selected geometry and
upstream entrance model; mixing junctions, bends, changing cross-section or wall
conditions are not newly modeled. Updating local properties or finite wall
and bulk temperatures is a held-property approximation. Exact partition
additivity is claimed only at the same source and wall state, and the original
mean-card allowance is not a validated segment-level model error bound.

Three integration regression sources cover independent whole-duct transfer,
heating/cooling and reversed entrance order; native coupled momentum, advection,
finite wall and reservoir energy; original gas film EOS/NASA calorics and
endpoint Mach refusal. Ownership, geometry, cancellation and crossing support
are exercised in those same tests. Source/API review was completed; Rust
parsing, compilation and native numerical execution were unavailable, so these
regression sources are not reported as executed.

### Source-integrated walls on developing duct sections

`developing::conjugate::PreparedGraphDevelopingConjugateLinks` binds the original
finite wall contact, identical wall caloric object and supported pure-phase
conductivity chart. Liquid and gas owners expose
`evaluate_developing_conjugate_links` and
`evaluate_quasi_steady_developing_conjugate_links`. Only existing zero-conductance
finite links are admitted; a prescribed reservoir has no finite material trace.
The normal-path area must exactly equal its original contact area, and the
section owner still requires the whole native wetted perimeter and length.

Ordinary and conjugate paths share geometry admission and the same original
circular/rectangular source evaluator. Geometry is admitted once per endpoint
before surface trials. Each actual surface uses the original bulk/mean-film
convention, signed forward/reverse entrance distances, initial/final/average
flow support, density and native gas Mach/work checks. The shared material
power-coordinate root balances supported k(T) integrals against actual section
film heat, retaining cumulative evaluations and the final source payload. One
Q/(Twall-Tbulk) coefficient transfers heat reciprocally; bare film hA remains
diagnostic. Exact contact/caloric identity and atomic conflict/cancellation
checks are shared with the other conjugate links.

The selected card retains its original inlet, uniform-wall and held-property
assumptions. Different cell/surface states do not establish exact partition
additivity or a resolved nonisothermal entrance history. The rectangular owner
still requires its original aspect ratio 0.5 and source Pr 0.72; no clipping or
substituted coefficient is introduced. Internal phase fronts, axial wall
spreading and new section-level model-error bounds are not supplied.

Three focused integration specifications cover independent Hausen/affine-k
surfaces, both heat directions and equilibrium, reversed entrances, original
rectangular table values, actual film-source temperature, native reciprocal
residuals, quasi-steady parity and typed physical/budget refusals. Native graph
tests were syntax/API reviewed and remain uncompiled/unexecuted with missing
constellation siblings. A separate bounded build passed 18 actual material and
rectangular arithmetic tests, including one coupled affine-k/section-heat
oracle; it does not compile the full source-card or native graph stack.

## Original gas correlations on planar moving-liquid heat interfaces

`graph::thermal::medium::convection::PlanarFreeSurfaceConvection` implements
fs-lbm's call-local surface coefficient provider. It consumes the original
source gas already resolved by the moving-liquid owner at the actual trial
mean-film temperature. Forced plate Reynolds/Prandtl groups and `h = Nu k/L`
use that same source state. Natural plates additionally use its original
volumetric expansion and the full actual interface/bulk temperature difference.
The original correlation evidence is moved into the accepted surface sample;
no surrogate gas, constant-h default or extra heat integrator is supplied.

The declared and actual liquid-to-gas normals must be finite unit vectors that
agree in world coordinates, including planar Cartesian, height-graph and PLIC
interfaces. Unit length, direction agreement and tangency admit only
`64 * f64::EPSILON` normalization roundoff, with no empirical angular allowance.
Forced flow is a prescribed gas velocity RELATIVE to the liquid surface, must be
tangent, and explicitly selects the existing laminar or turbulent average plate
card with its original regime limits. Natural gravity must describe a vertical
or horizontal plate; the outward normal and signed temperature difference select
the existing stable/unstable horizontal branch. Geometry, source and physical
driving inputs are admitted even at exact equilibrium, where natural heat is
zero without inventing a finite-Rayleigh coefficient.

`AidingVerticalPlate` admits laminar aiding mixed convection on that same
moving-liquid patch. Physical plate-up, actual gas velocity relative to the
liquid, gravity and a caller-admitted density-departure bound are explicit.
Plate-up and relative flow must be tangent to the actual face. The shared
material frame and group admission retain signed axial flow, source-owned
expansion and the full interface/bulk contrast; oblique gravity, crossflow and
opposing buoyancy are refused. The original Manning–Qureshi card supplies Nu,
its pure forced/natural limits, evidence and laminar regime bounds. Exact
resting equilibrium has zero heat without an invented Re=Ra=0 coefficient.
Mean-film sampling is this adapter's frozen-property approximation, not a
source-mandated rule; the explicit density bound does not establish weak
variation of viscosity, conductivity or other properties.

Characteristic length means the whole physical streamwise plate length, vertical
height, or horizontal area/perimeter. It is never reset to a mesh-face width.
Applying the mean isothermal-plate correlation to surface patches assumes an
approximately isothermal declared physical plate and locally held film
properties. Relative speed and geometry remain held during the original
thermal/native-flow split. Curved, enclosed and evaporating surfaces require
other physical models. Natural convection inclined to gravity and mixed
opposing/crossflow cases still refuse; no gas velocity or gravity component is
projected into an admitted orientation. The original liquid
half-cell, radiation, finite-gas reciprocal balance, thermal interval admission
and rollback remain owned by fs-lbm. The production adapter's unit tests compare
with the original correlation and independent plate equation, source-pressure
response, gravity orientation, equilibrium admission and explicit refusals.
World-rotation regressions cover both forced regimes, natural vertical and
horizontal branches in both heat directions, aiding mixed flow and its pure
limits, plus non-tangent flow and small but physical orientation mismatches.

The same single-plate provider implements fs-lbm's
`FreeSurfaceBodyGasConvection` for finite stationary solid faces. Both providers
delegate to one physical plate evaluator using the actual outward normal,
solved surface temperature and original bulk/film gas states. The body adapter
retains signed outer-wall geometry and never fabricates a liquid patch. Solid
conductivity, surface balance, finite gas storage and time integration remain
with the original fs-lbm/material owners. The same forced, natural and aiding
mixed cards, source support, physical length and orientation requirements apply.
Two focused body-provider test functions cover an independent forced-plate
equation, pressure dependence, both heat directions, natural/mixed evidence,
equilibrium and source/orientation refusal. They were syntax/API reviewed;
native compilation and execution remain pending constellation dependencies.

`PlanarFreeSurfaceConvection` also implements the stationary material/phase-flow
`FiniteGasConvection` trait. Its explicitly authored contact normal, solved
surface and original bulk/film gas pass directly to the same plate evaluator,
without fabricating a moving-liquid or finite-solid patch. Forced, natural and
aiding mixed laws retain original source, regime and orientation support. Two
provider specifications cover an independent forced-plate equation and pressure
scaling, actual film temperature, natural/mixed evidence, normal reversal,
opposing flow and equilibrium. Syntax/API review passed; native execution
remains pending an admitted constellation build.

The single-plate provider explicitly admits uniform boundary evolution: its
authored plate frame, dimensions and flow do not depend on old patch indices.
Coefficients still depend on the actual trial temperatures and original film
source. Newly exposed native patches must pass every material and actual-normal
check against the same declared plate; the capability does not admit curvature,
a different local flow or an unsupported natural-convection inclination.

### Streamwise external-plate sections

`evaluate_flat_plate_section_from_material` now supplies a section-average
coefficient from the original fluid state and the two physical distances from
one leading edge. It retains that source, its single transport snapshot, both
original cumulative correlation evaluations and the section coefficient.
The same coefficient lowers through the existing evidence-paired Robin
boundary into solid conduction. Laminar and leading-edge-corrected turbulent
cards keep their own domains at both positive endpoints; no transition blend
or new entrance at a cell boundary is introduced. A section starting at zero
retains the original whole-length mean exactly.

`CellConvectionLaw::ForcedPlateSection` and
`SourceCellConvectionLaw::ForcedPlateSection` carry the same section coefficient
into original finite-cell boundaries and reciprocal finite-fluid links.
`start` and the binding's `characteristic_length` are upstream and downstream
distances from one physical leading edge. Source bindings require mean-film
properties at the actual wall/fluid temperatures, and invalid geometry refuses
before source callbacks. Existing observers and heat assemblies apply the
coefficient to native exposed/contact area and include its conductance in the
original timestep bound. `CellConvectionEvaluation::section` retains the
upstream mean and physical interval; `nusselt` retains the original downstream
mean, while `effective_nusselt()` returns the actual strip mean. Positive-start
coefficients preserve the shared section owner's unbounded local-model
discrepancy. Regression sources exercise unequal subdivision, heating/cooling,
film refresh, native timestep bounds and closed two-store energy balances.

`PlanarFreeSurfaceSections` applies those coefficients to the original
liquid/gas heat solve by O(1) lookup in authored surface-contact order.
Every patch retains its liquid/gas identities, area, normal distances and
world-frame normal. Explicit downstream directions must agree with the
actual relative flow; reversed flow needs distances from the opposite leading
edge. Missing or stale patch geometry and unsupported upstream regimes refuse.
Film properties are re-evaluated by the original material owner at the actual
trial interface/bulk mean, and the native liquid resistance, radiation,
finite-gas capacity, reciprocal energy transfer and rollback remain active.
This spatially varying provider deliberately does not admit uniform boundary
evolution: newly exposed patches require their own physical section binding.

For a common isothermal wall, held fluid state and relative velocity, the
area-weighted strip coefficients recover the original whole-plate heat while
resolving its downstream variation. Locally re-querying film properties is an
explicit approximation to that constant-property boundary layer. Arbitrary
non-isothermal boundary-layer history, transverse variation, curvature and
gas momentum remain outside the relation. The source mean-card discrepancy
is not a local error estimate: positive-start section evidence carries an
explicit unbounded model discrepancy. Numerical checks cover unequal
subdivision, source/dimensional consistency, Robin lowering, local heat
variation, rotation, geometry changes and flow reversal; they do not confer
experimental validation.

### Plate sections on native hydraulic and thermal graph endpoints

`GraphConvectionLaw::ForcedPlateSection` connects the original external-plate
strip coefficient to source-bound liquid and gas graphs. Native solved `Q/A`
supplies the free-stream speed; `start` and `characteristic_length` retain the
two distances from the physical leading edge. `PlateFlowDirection` explicitly
declares whether downstream follows or opposes the branch path. Initial,
average and final flow must agree with that direction. An interval reversal
refuses before material queries; subsequent reverse flow requires rebinding
both direction and distances from the opposite edge.

Finite wall links and prescribed wall reservoirs use the same original
mean-film source and section evaluation. Both positive spatial endpoints must
pass the selected laminar or turbulent card at both temporal endpoint speeds;
a supported average cannot hide an unsupported upstream Reynolds number.
The native implicit graph refreshes these coefficients with its actual trial
temperatures and solved flow. Finite wall exchange cancels internally through
the original reciprocal link; only prescribed wall heat is external. No new
correlation, thermal storage, mass flow or source-power deposit is introduced.

Focused liquid/gas `source_graph::plate_sections` regressions cover held-state
strip additivity, explicit reverse-flow rebinding, original gas film/EOS and
caloric identity, interval/domain refusal, and a joint liquid momentum/enthalpy
step with viscosity and film-coefficient feedback and native energy closure.
The external plate's attached, unconfined, zero-pressure-gradient assumptions
remain caller obligations. Local temperature updates remain the existing
held-property approximation; they do not resolve boundary-layer history or
establish experimental accuracy. Native test execution is reported separately.

### Forced spheres on native graph endpoints

`GraphConvectionLaw::ForcedSphere` connects the existing Whitaker evaluator to original liquid and gas graph storage, using native solved speed and the whole sphere diameter. Bulk temperature, density, viscosity, Cp and conductivity supply the primary correlation state. The same pinned fluid source is additionally resolved at the actual wall temperature solely for its viscosity correction; source pressure, composition and identity remain unchanged. Finite-wall and prescribed-reservoir reports retain both material owners. The surface property query has no separate stored mass or energy, so only the primary bulk state receives the native reference-density admission.

The original sphere area, Reynolds, Prandtl and bulk/surface viscosity-ratio bounds remain active. Initial and final speeds use the retained two-state snapshot; a supported average cannot hide an unsupported endpoint or the zero-speed crossing inside a reversal. Equal already-admitted speeds avoid repeated correlation/evidence construction. Mean-film substitution, a unit viscosity-ratio default, extrapolation through stagnation, confinement and natural/mixed-convection corrections are not introduced.

The existing implicit native graph refreshes flow, bulk and wall properties at each actual thermal candidate. Original reciprocal links debit finite walls once; prescribed reservoirs enter external heat once. Whitaker's conductive contribution stays inside its single coefficient. Focused `source_graph::spheres` regressions cover nonunit liquid/gas viscosity ratios, source and interval refusals, gas surface-property versus storage-density separation, and coupled momentum/enthalpy evolution with finite-wall and total energy checks. The isolated stationary isothermal sphere and uniform external-flow assumptions remain; numerical execution and experimental validation are reported separately.

### Source-integrated solid resistance on native forced-duct links

`PreparedGraphConjugateLinks` composes the original liquid/gas graph convection
admission with `ResolvedConductionLayer` and
`ReservoirLayerGeometry::evaluate_material_to_bath`. One homogeneous,
source-free solid path runs from the existing wall cell's temperature trace to
its actual fluid-facing surface. Planar or radial geometry remains explicit;
its contacted area and original contact coverage must match exactly. The
prepared conductivity chart retains the original source line and identical
wall caloric object; no conductivity copy, additional solid storage, or
separate fluid temperature is introduced.

The original `new`/`resolve` fast path admits bulk-property circular/rectangular
fully developed laminar CWT/CHF, Dittus-Boelter, and Gnielinski laws only. The existing native
graph owner still checks the actual initial/final/average flow and
source/correlation limits. Film-property, sphere, plate, natural, mixed, and
developing laws refuse through that fast path. Surface-dependent native laws
use the explicit closure below; developing section bindings retain their
separate owner. The caller retains the selected correlation's developed-flow
and wall-uniformity assumptions.

At each nonlinear endpoint the existing material owner solves the surface
temperature and integrated k(T) heat. The adapter lowers its accepted signed
heat to the single original reciprocal link using Q/(Twall-Tfluid), with a
separate equilibrium limit and explicit representation audit. Only the
combined conductance map is applicable; bare film samples remain diagnostic
evidence. Reports retain the actual surface, original fluid source, shared
solid source line, and layer residual. Changed wall/caloric/contact ownership,
stale callback evidence, source-domain failure, competing conductance
ownership, and cancellation refuse without advancing either native history.

The focused `liquid_pressure::source_graph::conjugate` cases cover
constant/variable-k surface continuity, heating/cooling/equilibrium, native
reciprocal residuals, refusal cases, and coupled
momentum/advection/material-wall endpoint feedback. Ten tests passed in a
bounded build of the original material layer/reservoir kernel, including two
fixture-specific flux/domain checks. The new graph adapter and its integration
tests were syntax/API reviewed but not compiled in the full Cargo graph.

### Original source films at conducting-wall surfaces

`PreparedGraphConjugateLinks::new_with_surface_feedback` binds the same exact
wall caloric object, original conductivity source and contact geometry. Liquid
and gas graph owners expose `evaluate_conjugate_links` and
`evaluate_quasi_steady_conjugate_links` to solve the shared material integral
against their original film at each actual surface trial. This admits existing
film-property forced ducts, external plate sections, Whitaker spheres, natural
and aiding mixed laws without replacing any correlation or adding another
thermal cell. Mean-film properties use the solved surface and actual bulk
temperature; Whitaker retains its primary bulk tuple and separately resolves
only its original surface-viscosity state. Original geometry, source support,
signed native flow envelope, density, local-work and gas Mach limits remain
authoritative.

The shared material kernel returns the actual accepted film power, source
payload and solid/film residual. The graph publishes one Q/(Twall-Tfluid)
coefficient, preserving the explicit equilibrium limit; bare film coefficients
remain diagnostic and cannot be applied as another heat owner. Bracket/source
failure, exhausted power-residual budget, cancellation and competing overrides
refuse atomically. This is a quasi-steady homogeneous normal path with the
selected original correlation's uniform-wall and external-flow assumptions;
it does not resolve wall spreading or boundary-layer history. Developing-duct
sections use their explicit conjugate owner above.

Focused `source_graph::conjugate::surface` tests compare independent film
equations for a variable-k wall, duct/sphere/plate source states, natural/mixed
heating and cooling, equilibrium, native inertial/quasi-steady consistency,
reciprocal heat and refusal behavior. The existing full coupled
momentum/advection test now exercises the solved surface film. These graph
tests were syntax/API reviewed; full Cargo execution remains blocked by absent
constellation checkouts. Shared-kernel execution is reported separately.

## Reversible total-pressure hydraulic nodes and physical static tractions

`graph::inertial::total_pressure::TotalPressureGraph` is an explicitly selected,
fixed-geometry incompressible owner. `FixedTotalPressure` is distinct from the
original `FixedPressure`: it supplies `P = p + rho (Q/A)^2 / 2`, in the original
source graph's pressure datum. A stagnant reservoir may prescribe this head.
Finite-speed static boundary pressure is not silently converted.

The original descriptor retains each branch's impulse `rho L Q/A`, quadratic
kinetic storage, exact zero-storage incidence constraints, held drag and
hydrostatic reduction. Each branch has uniform area, so the same dynamic head
is subtracted at both of its ends and cancels from its impulse equation.
No second momentum solver, pressure penalty, invented capacitance or steady
flow root is introduced. The original default static-pressure graph is unchanged.

Every incident port has its own static pressure. For the native affine
trajectory, mean `Q^2` sets interval-average traction, while mean `Q^3` sets
pressure/kinetic work. The existing stable formulas are shared:
`mean(Q^2) = qbar^2 + dq_half^2/3` and
`mean(Q^3) = qbar (qbar^2 + dq_half^2)`.
Multiplying average static pressure by average flow would give the wrong work
during acceleration or reversal; the new path does not use that shortcut.

`BranchStreamWork` reports both physical port pressures and signed static,
kinetic and total-pressure work. `NodeStreamWork` retains each complete
incident set. At a reversible zero-storage junction, its pressure work supplies
the outgoing kinetic-energy transport, including unequal areas and reversal.
The full graph report independently closes original kinetic storage and loss
against actual external static-pressure work, kinetic-energy inflow, explicit
drives and full physical gravity. Local node residuals and pressure/kinetic
identity defects are kept distinct from the native descriptor's own ledger.
Both those residuals and the original native impulse/work gates must pass
before acceptance.

The original pipe, junction and assembly load consumers use these actual
static tractions. Internal pressure and advective cuts cancel before the
whole-assembly sum; fluid acceleration, full bend angular momentum and weight
remain exactly once. `accept_with_assembly_loads` publishes the same accepted
native trajectory and its usable load package together. Explicitly enclosed
pressure-datum nodes still require complete endpoint flow continuity and their
own local work acceptance; external streams cannot disappear by declaration.

Material-bound pressure admission checks both static-pressure extrema over the
affine interval, including its interior maximum when flow crosses zero.
Original reference-density, source velocity/pressure support, topology,
gravity, nonlinear closure and neglected thermal-work limits remain in force.
Rechecking refreshed coefficients audits the same trajectory. Cancellation,
unsupported source state, incomplete geometry or failed work acceptance leaves
the initial owner untouched. The bounded native state cap also bounds the
additional ordered stream reporting; no field resampling occurs.

The total-pressure owner keeps its raw native report private, so its total
multipliers cannot be passed to existing static-pressure thermal/source
adapters under the wrong chart. This increment supplies a hydraulic momentum
and stream-work closure; it does not add finite-speed static-pressure boundary
conversion, turbulent junction mixing, a dissipation-to-heat credit, acoustic
storage, variable-density or cavitating flow, swirl, nonuniform section
profiles, or moving-pipe FSI. Those require their own physical owners.

## Unsteady circular laminar flow with viscous memory

`graph::inertial::laminar::LaminarPipe` binds the original `FluidTransport`
owner and one frozen source snapshot to the existing circular resistance and
positive `fs-phs::ZkDiffusiveBasis`. Its native PHS contains mean-flow inertia
and all admitted viscous-memory states, with separate profile storage and
irreversible wall loss. Steady initialization recovers Poiseuille resistance
and the parabolic profile's 4/3 kinetic-energy factor. Memory pressure work is
never automatically credited as heat.

`PreparedLaminarPipe` uses the original exact held-input trajectory, integrated
port flow and independent native heat integral. Preparation caches the bounded
dense work; successful warm steps and mean-flow relations allocate no memory.
Two pressure ports conserve volume. Native-system handoff supports original
PHS composition and transfers source/history/regime admission to that assembly.
A comparison envelope admits the entire standalone held interval's mean-flow
Reynolds support. Source, history-owner, regime, arithmetic, native-solver and
cancellation refusal leave accepted history unchanged.

`LaminarPipe::with_source` continues a different original material snapshot on
the same physical pipe without resetting mean flow, viscous profile memory,
kinetic energy or accepted time. It retains the original dimensionless diffusion
basis and pole ordering; no eigensolve is repeated. Actual source viscosity
rescales all viscous gains and rates together, leaving their storage weights
unchanged. Every proposal samples its source once and repeats the original
frequency-band complex-response and wall-loss admission. A newly prepared
propagator evaluates the continued history under that source's Reynolds support.
Old owners/preparations remain immutable. Failed source, band, regime or
arithmetic proposals leave the original history available for continuation.

Geometry and density are fixed; even a one-ULP density change refuses because
altering the finite pipe's mass/kinetic storage needs its own conservation owner.
Viscosity is spatially uniform and held on each accepted interval; the caller
owns the thermal or composition path and time refinement. This continuation
does not infer temperature evolution, radial property gradients or fluid mixing.
Focused G1/G3 regressions continue a developing profile through a source change,
compare with independently initialized equivalent Stokes evolution, recover the
viscous-time rescaling of unforced decay and heat, and exercise density, renewed
band and Reynolds refusals with untouched accepted history.

Assumptions are stationary rigid no-slip walls and axially developed
incompressible Newtonian flow. The original
viscous spectrum retains sampled frequency-band/error admission; arbitrary
startup has no inherited time-domain error certificate. Entrance effects,
turbulence onset, radial heat transport, cavitation, compressibility and
changing density require other owners. Gas use additionally needs
caller-qualified low-Mach/frozen-density support. Native regressions exercise
independent radial Stokes startup/reversal and kinetic energy, harmonic phase
and power, finite-reservoir composition, and atomic history/regime admission.

### Continuous pulsatile pressure histories

`LaminarPipe::prepare_tones` drives that same material-bound pipe and viscous
history with simultaneous continuous pressure-difference tones, including DC.
The original `PreparedExactTones` owner integrates the full forced trajectory,
displaced volume, pressure work and irreversible heat, retaining cross-frequency
terms. Pressure is not held at a sampled value inside the interval. This removes
forcing-hold error without changing the physical wall spectrum or dropping
profile memory. Every oscillatory frequency must be inside the original sampled
wall-response band; DC retains the exact Poiseuille limit.

`PreparedLaminarTones::advance` advances the pipe clock and caller-owned source
phase coefficients atomically. Restarts retain both. The sum of pressure-tone
amplitudes supplies a conservative comparison envelope for the entire interval's
Reynolds admission, including reversals and peaks missed by endpoint samples.
It may refuse strongly cancelling waves; it is not a turbulence predictor or a
time-domain continuum error certificate. Original fixed geometry/density,
Newtonian, developed, incompressible and rigid no-slip assumptions still apply.

Preparation owns the explicit dense work/memory budget and reusable buffers.
Warm calls allocate nothing and use cached phase rotations. G1/G3 regressions
compare pulsatile startup with independent radial Stokes modes, recover the
held-pressure DC path, preserve flow/displacement/work/heat under interval
subdivision, and retain accepted flow and forcing phases on refusal.

`LaminarNetwork::prepare_tones` extends this same continuous forcing to the
complete constrained network. Every original reservoir and directed branch
pressure-head port supplies the full prepared frequency list in port-major
order. Each connected component's pressure gauge is removed separately for
each quadrature before matrix products. All original radial-memory arms remain;
internal junction continuity is built into the original coordinates.

The native trajectory independently integrates every original branch's positive
viscous source squares and every physical port's signed pressure work, including
cross-frequency effects. Reservoir work and branch-head work remain separately
available, with positive work entering the fluid. Individual reservoir-port
work uses the component's reference gauge; the total reservoir work is gauge
invariant. The adapter audits branch heat against native total heat, source work
against native supplied work, and endpoint and integrated reservoir/junction
continuity. It publishes history, clock, source phases and all accepted reports
atomically. Failed admission or cancellation preserves them for retry.

Preparation accounts for complete channel construction and all reusable buffers
within the caller's explicit dense byte/work/dimension budget. Warm advancement
allocates nothing and performs no new exponential or factorization. Passive
contraction bounds the entire interval's Reynolds support; a separately retained
DC equilibrium sharpens the bound. The oscillatory forcing-norm envelope is
conservative and can refuse long intervals or strongly cancelling waveforms.
The caller may request a shorter interval; the owner never changes the clock or
silently removes modes. Developed, fixed-density Newtonian branch assumptions,
zero-storage junctions and original wall-spectrum validity limits still apply.
These ports prescribe pressure, not a pump curve or an automatic flow controller.

Focused regressions compare a one-branch network with the original continuous
pipe, compare a fork's full trajectory and every heat/work channel with an
independent RK4 integration, preserve history and all accounts under subdivision,
recover held-pressure DC behavior, preserve independent component gauges, and
retain accepted reports/history/phases when an interior pressure peak or bad
input is refused.

## Material-backed finite rectangular cooling-channel sections

`graph::thermal::medium::convection::duct` evaluates the original rectangular
CWT developing table over physical heated-entrance distances. One original
`FluidTransport` snapshot supplies Re, Pr and the conductivity in the section's
dimensional h. Hydraulic diameter, aspect ratio and both distances remain
explicit physical inputs. The existing source still admits exactly Pr=0.72
and aspect ratio 0.5; other fluids/geometries and unsupported positive upstream
positions are refused. Opposite flow requires authored opposite-entrance
distances. No cell boundary silently resets development.

The source-bearing section coefficient lowers to the existing conduction Robin
row with the actual bulk reference temperature. Uniform-state subdivision
preserves cumulative heat; varying local states remain a held-property
approximation, and positive sections have no validated local discrepancy band.
Two source/material integration regressions cover conservation and original
source-domain refusal; they require the unavailable constellation closure and
are not reported as executed here.

## Source-resolved developed sections in series

`graph::inertial::laminar::series::LaminarSeries` connects the original admitted
`LaminarPipe` sections with exactly one common incompressible flow. It retains
each source, radius, viscosity and every wall-memory arm. The duplicate mean
flows are eliminated analytically, with total inertance equal to the sum of the
original section inertances; no intermediate compliance, steady-flow root or
step-average continuity approximation is introduced. One energy-normalized mean
coordinate and all original profile coordinates carry the exact sum of the
section kinetic energies. Construction checks the complete 128-state bound and
exact common reference density before allocating dense matrices.

The original native exact propagator advances the entire held-pressure coupled
trajectory. Its port average drives displaced volume and work. Independent
`PreparedSourceLoss` caches integrate each section's original positive viscous
source squares on that same trajectory, preserving local heat deposition even
while the profiles develop or reverse. These integrals close against the native
total dissipation before state and caller-owned section-heat slots publish
together. They are never assigned by a steady resistance fraction, an endpoint
power estimate or a kinetic-energy deficit. Preparation owns all exponentials
and scratch; successful warm advances and mean-flow relations allocate nothing.

`with_sources` renews each original pipe's source, frequency-band and viscosity
admission while preserving the coupled mean flow, every local profile history,
kinetic energy and clock. Every density must remain exactly unchanged. The
cooperative mean/profile-flow comparison bounds Reynolds support for every
section throughout a held interval; each original section's ceiling applies.
A failed source, owner, regime, arithmetic, heat closure or cancellation leaves
accepted state and local heat output unchanged.

This is an ideal series connection of developed Stokes components. Different
radii do not infer area-change, entrance or convective momentum losses, and
local heat output does not itself evolve a thermal inventory. Thermal expansion,
radial viscosity variation, density changes, mixing and turbulent transitions
require their own owners. Tests check uniform-pipe subdivision through startup
and reversal, each section's heat and kinetic inventory, heterogeneous viscous
history against the original full common-flow descriptor, unchanged energy
through a nonsteady material update, steady local heat, complex impedance and
atomic physical-limit refusals. Runtime verification is reported separately.


### Particle-size-resolved soot opacity and direct-beam transport

RayleighSoot::particle_opacity and particle_opacity_from_mass_concentration
prepare a source-borrowing monodisperse sphere report with distinct absorption,
scattering and extinction in inverse metres. The original soot spectral
absorption is retained exactly; elastic Rayleigh scattering is added only to
primary-beam extinction. Source identity, condensed-particle density, radius,
vacuum wavelength and original population remain available from the report.
The host is nonabsorbing with refractive index one.

Admission requires the original loading domain and positive finite size and
wavelength, with x <= 0.1, |m| x <= 0.1 and x² |(m²-1)/(m²+2)| <= 0.1.
These are explicit asymptotic-use restrictions, not certified approximation
error bounds. Aggregation, particle-size distributions, dependent scattering
and the full Mie regime require separate models. The source must justify its
optical constants at the requested wavelength.

RayleighParticleOpacity::direct_beam returns unscattered transmission,
expm1-evaluated thin-path opacity, and separate primary-beam absorbed and
scattered fractions. Their energy partition and consecutive-path composition
are checked against analytic references. Scattered power remains radiation:
the returned absorbed fraction is first-interaction absorption, not total heat
after multiple scattering. Existing thermal LayerAbsorption and Planck means
receive no scattering coefficient. Preparation and beam queries allocate
nothing, and overflow of a finite positive optical-depth product has the
opaque limit. This is not a human-visibility or full radiative-transfer solver.

## Diffuse transmission through participating-enclosure boundaries

`radiation::sheets::SheetPair` binds two equal-area, equal-temperature faces
to one source-bearing `ReciprocalSheetOptics`. `evaluate_with_sheets` and
`evaluate_room_with_sheets` reuse the original reciprocal optical chords,
gas absorption/emission, workspace and LU solver. Gray material cards are
supplied only for unpaired opaque faces, in increasing view-factor order;
sheet sources exclusively own the paired optical coefficients.

Each outgoing sheet face contains its own emission and reflection plus
irradiation transmitted from the opposite face. The original exchange matrix
retains implicit self-return. Thermal face powers are actual absorbed minus
emitted power; transmitted energy is never deposited in the sheet. Gas and
opaque-boundary powers retain the original reciprocal difference operator.
The report preserves the original sheet and opaque sources, physical equation
residual and measured whole-enclosure energy closure, without redistribution.

A directed dependency check rejects closed lossless optical components before
factorization, including cycles formed by transmitting sheets. Pair overlap,
unequal physical areas or temperatures and unsupported optical temperatures
refuse. The no-sheet path retains the original opaque assembly and arithmetic.
This is an isothermal diffuse, same-exterior-index model; specular windows and
through-thickness temperature gradients require their own transfer model.

### Wavelength-selective reciprocal sheet boundaries

`SpectralOptics::new_with_sheet_bands` binds one sourced sheet partition to
each existing gas/surface wavelength band. Every band must retain the same
ordered face pairs. Its original full-axis surface absorptance must match the
paired sheet source exactly; conflicting sources refuse. Both temperature
domains and the enclosure's equal-area/isothermal pairing admission apply at
every evaluation.

Band-constant optics use the original integrated Planck emission differences
and one matrix factorization/solve per band per evaluation; a zero point-wavelength
sampling budget is valid for this exact-band path. Wavelength-varying
soot uses the original bounded adaptive quadrature and the same generalized
sheet operator at each sample. No duplicate spectral integrator or inferred
gray spectrum is introduced. Finite-window coverage remains explicit and
unrenormalized; nonabsorbing faces create no fictitious emission source.
Original sheet sources, gas sources, numerical work, spectral diagnostics and
actual material/gas heat remain in the returned report.

### Implicit original-inventory heating of transmitting sheets

`radiation::sheets::thermal` advances the original
`fs_thermal::enthalpy::DepositionSimulation` through its native
`BackwardEulerStep`. Both optical faces of each sheet must bind the SAME
active `CellId`. Their actual absorbed-minus-emitted powers are summed once
into that cell's endpoint source. No sheet mass or caloric law is duplicated.
Unpaired opaque faces may bind existing cells or explicit prescribed baths.

The candidate solves energy increments for every original active cell, also
retaining cells without an optical face, original conductive links, material
domains, authored sources and latent plateaus. Immutable old history is reused
for all nonlinear probes. A bounded dense Newton workspace admits at most 256
cells; source-domain refusals shrink finite-difference probes and line-search
steps under explicit evaluation/iteration limits. Local and total native joule
residuals, plus the actually applied cell/gas/bath energy balance, gate
acceptance. Cancellation or refusal leaves the input thermal inventory intact.

The returned material endpoint already includes the reported cell radiant
powers. Apply its original gas powers once and retain the explicit prescribed
bath debits when accepting the coupled owners. Gas temperatures, optical
geometry and sourced sheet coefficients are held over this interval; gas is
not advanced here. Native thermal clock/event admission remains in force:
interior events refuse and accepted endpoint events remain pending until
joint-owner acceptance. These are algebraic and energy checks, not a temporal
accuracy or through-thickness isothermal certificate.

## Non-isothermal absorbing material in participating enclosures

`radiation::slabs::SlabPair` binds two equal-area exterior ports to one
original prepared `AbsorbingSlab`. The gray and spectral slab entry points
accept one actual temperature per original material layer; each port matches
its adjacent first/last layer. Front and back temperatures may differ.
The original gas paths, view factors and factored boundary equations continue
to own enclosure exchange. The opaque path retains its original arithmetic.

Material heat is reported separately in stack-then-layer order. Exterior ports
have exactly zero surface thermal power: their incoming-minus-outgoing flux
must not be deposited alongside actual layer absorption minus emission.
The full energy gate sums opaque boundaries, gas and every material layer.
The same directed dependency admission rejects indeterminate lossless cycles.
The prepared material response retains angular memory inside each stack;
exterior coupling uses the enclosure's explicit diffuse approximation.
Sourced uniform or heterogeneous Fresnel slabs retain each face's reflectance in
the factored boundary equations and surface residual. Reflected power feeds
the enclosure, while internally trapped modes remain in the original layer
heat coefficients; neither is deposited a second time at the exterior ports.
The same gray, spectral and native thermal paths consume this response.
Distinct layer indices and asymmetric front/back reflection stay in the original
slab owner; the enclosure does not replace them with an averaged reflectance.
The boundary operator retains its original H=J+P/A form: the slab diagonal is
transmission plus that face's summed layer absorption, with reflected self-return
added once by the shared factorization. This positive 1-R representation also
preserves weak escape when reflection rounds to one.

`SpectralOptics::new_with_slab_bands` requires genuine band sources for
non-port faces and each physical material layer. Every band retains identical
ordered ports, layer count/order and physical thicknesses, while absorption
and source support may vary. Every emitting layer contributes its own Planck
spectrum and spectral-coverage gate, including layers hidden from both ports.
Constant band optics need one factorization/solve per band and no wavelength
samples; varying soot uses the existing adaptive wavelength integration.
Reported layer heat errors cover wavelength quadrature. Each prepared slab
separately retains its angular coefficient diagnostics; source and spatial
errors are not inferred from either estimate.

`radiation::slabs::thermal::candidate_slab_thermal` binds these ordered
layers to active original `CellId` values. Each layer in one stack has a
distinct cell; separate optical stacks may share a cell and their powers are
summed without duplicating its mass. Both ports bind their adjacent layer
cells. Other faces retain stored-cell or prescribed-bath ownership.

The slab adapter shares the existing sheet adapter's bounded native
backward-Euler/Newton driver. All original cells, conductive links, authored
sources, material domains, latent storage and event/clock rules remain active.
Each immutable-history trial evaluates actual endpoint layer temperatures.
Gray callbacks query only non-port faces; genuine spectral inputs skip them.
Only actual layer and opaque-face heat enters native storage. Local and total
native joule residuals and the applied material/gas/bath balance gate acceptance.

`candidate_slab_thermal_with_feedback` additionally resolves sourced layer
absorption and real refractive index inside every native endpoint residual,
Jacobian and globalization probe. `SlabLayerState` carries the original
reference layer, original slab/layer/band axis and finite wavelength endpoints,
plus the actual native cell, material identity, specific enthalpy and temperature.
Thus latent-plateau phase response need not be misrepresented as a temperature
change. The deterministic caller-owned constitutive law must not advance its
own accepted history or infer spectral data from a total-gray coefficient.

Changing optics rebuilds the original incoherent angular/Fresnel response,
including heterogeneous refraction, total internal reflection, n² internal
emission and repeated returns. Identical layer responses reuse the immutable
coefficients. Each preparation retains its original quadrature budget and the
whole solve retains the shared thermal trial ceiling. Thickness, ports and
ordered cells remain fixed; changing geometry requires its physical owner.
The returned layer's source domain must cover the actual queried temperature.
Spectral feedback retains the original finite bands, opaque surface sources,
gas/soot sources and omitted-emission gates. Reports retain the optical sources
at the accepted thermal endpoint, never a stale initial-temperature response.

The returned candidate already contains its reported material heating. Apply
the returned gas/bath exchanges once when accepting their owners. Gas
temperatures and geometry remain prescribed; the fixed-optics entry point also
keeps its original coefficients fixed. Neither adapter advances gas or certifies
timestep, spatial resolution, scattering, coherence, phase kinetics or material
accuracy. A source-domain, geometry, budget or cancellation refusal does not
advance the input thermal state. Focused regressions compare sensible and
latent feedback with independent Fresnel ray sums and scalar energy roots,
and finite-band feedback with independently integrated band heat.

## Developing turbulent tube material coupling

`evaluate_turbulent_tube_section_from_material` carries the original material
owner, one bulk transport snapshot, the actual circular diameter and both
distances from the physical heated entrance into the shared Gnielinski
developing-section evaluator. The same density, viscosity, heat capacity and
conductivity determine Re, Pr and the dimensional coefficient. Native Robin
lowering retains the section coefficient and its evidence without another
material query. Finite section transfer is additive under fixed properties
and thermal driving; a new mesh cell does not restart the entrance. Reversed
flow requires rebinding the physical distances.

Source-backed `Forced` cell and link bindings also admit the
`GnielinskiDeveloping` mean with explicit heated L/D and circular geometry.
They require `Reservoir` bulk properties, exactly matching the correlation's
constant-property K=1 specialization. Film substitution and rectangular
geometry are refused before resolving material states. Wall viscosity is not
queried or silently used as a correction. Bulk state changes can be coupled to
the existing conservative finite-cell thermal step, while strong spatial or
cross-section property variation remains outside this held-property model.
The tube tests include independent corrected-formula references, physical
entrance subdivision, material query counts, native Robin lowering and closed
finite-body heat exchange in both link orientations.

## Turbulent entrance transfer on native coolant graphs

`DevelopingDuct::TurbulentCircular` selects the corrected constant-property
Gnielinski developing relation through the existing `GraphDevelopingSection`
owner. Native pipe diameter, flow area, section length and full wetted-wall
coverage must agree. The actual solved volume flow and the retained original
bulk transport determine Re, Pr and conductivity; `Reservoir` property
sampling is required by the explicit K=1 specialization.

Forward and reverse flow use their respective authored physical entrance
distances. Both endpoints and the interval-average native flow must remain
inside the original correlation domain; an interval reversal crossing zero
cannot hide behind an admissible average. Finite-wall links, prescribed
reservoirs and source-integrated conducting wall layers reuse the same section
evaluation, source/density admission and original thermal storage. Section
conductances are applied once through the native link/boundary maps; diagnostic
heat is not an extra source. Fixed-property partition additivity remains
conditional on one isothermal passage and common bulk state. The
`liquid_pressure/source_graph/developing` regression checks solved flow,
upstream/downstream transfer, reversed entrance positions, bulk convention and
the unsupported endpoint regime with independently evaluated formula values.

## Physical sections in native implicit and finite-cell heat exchange

`SourceCellConvectionLaw::ForcedTubeSection` and its borrowed-material
counterpart retain the actual circular diameter and both distances from the
heated inlet. They use Reservoir bulk properties for the shared developing
Gnielinski K=1 model. Native exposed-cell boundaries, reciprocal finite-fluid
links and shared heat/circuit consumers retain one thermal exchange and the
source at the current trial temperature. Tube patch area cannot exceed its
physical wetted interval. Returned `tube_section` metadata and
`effective_nusselt()` distinguish the actual section coefficient from its
retained original downstream mean.

`ImplicitLiquidConvection` also admits finite plate strips and tube sections.
Prepared plate-film property curves provide the full temperature chain through
the actual section coefficient's Re, Pr and conductivity partials, including
the half-factor from mean-film temperature and density conversion for kinematic
viscosity. Independent fixed-group boundaries do not remove available moving
group derivatives; a required missing partial or genuine source kink refuses.
Turbulent downstream strips differentiate their cancellation-corrected section,
not the whole-plate mean. Bulk K=1 tube coefficients have zero surface-
temperature derivative. Both retain upstream/downstream source evidence and
their actual section Nu, while native surface elimination, conduction,
radiation, body storage and time integration remain the existing owners.
Focused regressions cover material queries, signed reciprocal heat, entrance
position, analytic film chains and native body energy/entropy closure.

### Source-resolved liquid wall correction on the native surface

`ImplicitLiquidConvection::with_gnielinski_wall_correction` explicitly selects
the shared liquid `(Pr_bulk/Pr_wall)^0.11` correction for an existing developing
turbulent circular-tube section. Bulk Re, Pr and conductivity stay with their
original reservoir state. Wall Cp, conductivity and viscosity come from the
same prepared source at the actual eliminated surface temperature, including
the density factor for kinematic viscosity. Their complete analytic wall-Pr
chain enters the original nonlinear boundary; it is not evaluated at the stored
bulk-solid temperature. `surface_transport` retains that independent property
sample. The corrected section Nu/coefficient remain distinct from its retained
original constant-property cumulative means.

The prescribed-mass-flow wrapper retains this selection and the actual corrected
section Reynolds partial, alongside its separate wall-temperature tangent.
Source support, genuine source kinks and the correlation's property-ratio domain
remain active. Unselected laws keep K=1; this option adds no wall storage,
energy source, turbulent closure or second surface solver. Focused regressions
check both viscosity bases, independent corrected-section algebra and tangents,
held-rate flow response, source refusal, and original conductive/radiative
surface, backward-Euler energy and bath-entropy balances.

### Explicit air wall correction on the native surface

`ImplicitGasConvection::with_gnielinski_air_wall_correction` selects the shared
Gnielinski air factor `(T_bulk/T_wall)^0.45` on a developing turbulent circular
tube. Eligibility requires the retained original dry-air source; humid air,
named mixtures and thermochemical gases are not inferred to satisfy the air
model. Both bulk and actual eliminated wall temperatures retain declared
source support and use absolute kelvin. Original bulk transport still supplies
Re, Pr and conductivity. The native surface tangent includes `-0.45 h/T_wall`;
the coefficient also retains the independent `T_bulk` partial at fixed transport
groups, without claiming a total bath-temperature Jacobian.

The shared response records explicit turbulent correction inputs separately
from original K=1 cumulative means. Prescribed mass flow carries either the
liquid or air selection into its actual corrected section Reynolds partial.
No wall transport curve, storage or additional heat solve is introduced for
the air temperature factor. Focused G1/G3 regressions cover correction values,
source/geometry refusal, both analytic temperature axes, actual bulk-density
flow response and the existing native conduction/storage/bath-entropy balance.

## Mean laminar tube wall viscosity in source and native heat exchange

Selecting `SourceCellConvectionLaw::Forced` with `SiederTateLaminar` now retains the original mean laminar entrance correlation through both material callbacks and the native implicit exterior boundary. `Reservoir` denotes the authored arithmetic-mean bulk temperature for this model; the explicit L/D is the whole physical heated tube length. Bulk Re, Pr and conductivity stay at that source state, and only the same fluid's dynamic viscosity is queried at the actual wall temperature. The source's open viscosity/Pr intervals, laminar Reynolds support and coupled thermal-entrance gate remain active.

Source callbacks lower through `CellConvectionLaw::ForcedLaminarTube` and the existing material mean evaluator. Positive finite-cell area must equal the complete wetted area pi D L; a remaining patch cannot inherit the whole-tube mean. The pointwise native interface documents the same full-area requirement for its downstream boundary geometry. Neither path invents an axial profile or infers inlet/outlet temperatures from a local fluid state.

Liquid and gas native boundaries reuse their original analytic wall-viscosity queries. The coefficient tangent is `-0.14 h (dmu_wall/dTs)/mu_wall`, including the original density factor for kinematic liquid viscosity; the existing heat tangent adds `h + (Ts-Tmean) dh/dTs`. Prescribed mass flow retains the original laminar Reynolds partial. Source kinks or unavailable active partials refuse. Native conduction, radiation, storage and heat/entropy accounting remain in the existing solver.

Focused regressions cover dynamic/kinematic liquids, original air and CO2 sources, wall and mass-flow tangents, retained source coordinates, full-area admission, model/source refusals, and native surface/storage balance against an independent whole-tube law. These check implementation and conservation; they do not supply experimental validation or an error enclosure.

## Developing tube sections on finite gas heat owners

`FilmGasHeatLaw::BulkForcedTubeSection(ForcedTurbulentTubeSectionGasSurface)`
carries physical heated-entrance distances and the actual circular diameter
through the existing finite-gas/body, piston, fixed-volume wall-feedback and
mixed heat-network callbacks. Patch area cannot exceed pi D (end-start), and
flow reversal requires rebinding the physical distances. Each active callback
resolves one source-pinned bulk property state at the actual trial gas
temperature, pressure and composition, retaining the original ideal or explicit
real-gas branch. The shared developing Gnielinski evaluator uses K=1 without a
film-temperature or wall-viscosity correction.

`GasConvectionEvaluation::tube_section` retains physical geometry and upstream
mean evidence; `effective_nusselt()` returns the actual section Nu while
`nusselt` remains the original downstream cumulative mean. Existing finite
heat owners retain one signed reciprocal exchange, body enthalpy, pressure
work, area/Biot admission and atomic publication. Section additivity requires
common bulk properties and thermal driving; local refresh does not resolve
spatially nonuniform entrance history. Original correlation domains and the
modest-property-variation assumption remain active. The ordinary
`BulkForced` path also admits the developing mean with explicit circular
geometry and heated length. Regressions cover subdivision, source identity,
piston/body and fixed-volume energy, latent heat, and mixed heat-network
exchange.

`evaluate_laminar_tube_section_from_material` and
`FilmGasHeatLaw::BulkForcedLaminarTubeSection(ForcedTubeSectionGasSurface)`
carry the original Hausen thermal entrance into the same material Robin and
finite-gas heat owners. The circular velocity profile must already be developed
at the heated inlet. Both positive endpoint source domains remain mandatory;
neither a thermal patch nor flow reversal silently creates a new entrance.
One original bulk snapshot supplies Re, Pr and conductivity. Circular patch
area cannot exceed pi D (end-start). `TubeSectionDetails` retains actual section
Nu and the upstream mean; the previous turbulent metadata and geometry names
remain compatible aliases. The downstream mean is never relabeled as a local
correlation, and positive-start coefficient evidence retains unbounded local
model discrepancy.

`FilmGasHeatLaw::BulkForcedRectangularDuctSection` uses the existing
material-backed simultaneously developing rectangular CWT section. The source
still admits exactly aspect ratio 0.5 and actual Pr=0.72; no property value is
adjusted to fit that slice. A section must cover all four isothermal walls:
area equals Dh(1+alpha)^2/alpha times (end-start), within floating-point
roundoff. The source's perimeter-average coefficient cannot be applied to an
isolated wall or partially heated perimeter. `rectangular_section` preserves
hydraulic diameter, aspect ratio, physical distances, upstream mean and actual
Nu, including the original table/engineering-bridge distinction.

Both new film laws re-query the source at actual bulk temperature, pressure
and composition on every thermal trial, retaining original source and branch
admission. Zero-area isolation bypasses property/correlation evaluation.
Piston, fixed-volume and existing heat-network callbacks retain their original
single reciprocal heat transfer and finite storage. No second heat integrator
or frozen-temperature reservoir is introduced. Focused laminar regressions
exercise independent formula/table subdivision, source and geometry refusal,
piston heat in both directions and through melting, and both shapes on the
fixed-volume body/gas energy residual with external heat budgets.

The finite-cell adapters admit explicit `ForcedLaminarTubeSection` and
`ForcedRectangularDuctSection` laws on both `CellConvectionLaw` and
`SourceCellConvectionLaw`. Reservoir boundaries, reciprocal solid/coolant links
and the existing complete heat/circuit owner share the same material evaluation.
Each source-backed residual resolves the original fluid at its actual bulk
temperature; neither law samples film temperature or applies a wall-viscosity
correction. Both positive endpoint domains and the physical entrance remain
mandatory. A laminar query outside the original source regime refuses instead
of switching to a turbulent card. Fully developed inlet velocity is still
required by Hausen; the rectangular card retains its simultaneously developing
assumption and exact aspect ratio/Pr support.

Circular exposure cannot exceed pi D (end-start). Every active rectangular
binding must own the complete four-wall wetted area, within roundoff; partial
wall heating and accidental use of one face's area refuse. Native contact area
and exposed area remain authoritative. `CellConvectionEvaluation` retains
`tube_section` or `rectangular_section` geometry and the original cumulative
mean; `effective_nusselt()` returns the actual section coefficient. Positive
starts retain unbounded local model discrepancy. Native finite-cell storage,
reciprocal conductance, external heat accounting and trial acceptance are
unchanged. Regression tests close both caloric stores against independent
Hausen and rectangular table-intercept calculations, including reversed link
endpoint ordering and explicit source/area refusal.

`ImplicitLiquidConvection` also admits bulk-held Hausen circular sections on
the existing conducting surface and backward-Euler owner. Because every source
coordinate is independent of surface temperature, its coefficient has exactly
zero surface-temperature derivative; its heat tangent remains h, including at
thermal equilibrium. The native half-cell resistance, radiation, source caloric
curve and entropy ledger still participate in the same solve. The pointwise
exterior interface explicitly refuses rectangular four-wall sections because
it cannot validate their integrated wetted geometry. No empirical local law,
surface store or replacement heat integrator is introduced.

### Whole-tube laminar selection on finite gas heat owners

`ForcedTubeGasSurface::laminar()` explicitly selects the original Sieder-Tate laminar whole-tube mean for the existing finite piston/body and fixed-volume gas/body heat owners. A positive exchange covers the full physical heated area pi D L, permitting only geometric floating-point roundoff; zero area remains exact insulation. The original homogeneous bulk gas represents the caller's arithmetic-mean inlet/outlet bulk state, while only the same source's actual wall viscosity enters the correction. The whole heated length and source/laminar/thermal-entry limits remain active at every nonlinear trial. No downstream section profile, axial flow solution or automatic turbulent fallback is inferred. Both original energy stores, piston work where present, body phase change, wall feedback and entropy accounting stay with their existing native solvers.

The explicit laminar wrapper and default turbulent interface share the same material resolution and coupling routines. Focused regressions cover an independent whole-tube mean and source/geometry admission, reciprocal piston/body heat including latent enthalpy, and fixed-volume wall feedback with energy and entropy balance. These are implementation and conservation checks, not empirical gas-validation evidence.

## Pump drive input and coolant heating

Fixed-speed liquid cooling pumps may bind an independently sourced `PumpDriveCurve`
through `PumpHeadCurve::with_drive_curve`. Input samples specify pump-inlet flow,
shaft power and electrical power at the head curve's same fixed speed, with
explicit inlet-density and kinematic-viscosity support. Bounded piecewise-linear
interpolation admits only actual trials satisfying
`electrical >= shaft >= rho_in * g * H(Q_in) * Q_in`. This check is also applied
between head samples, where hydraulic output is quadratic in flow. Nonzero
shutoff consumption is retained; no constant efficiency, affinity scaling or
unprovided motor data is inferred.

`PumpLossRouting` separately declares the fraction of pump loss
`shaft - hydraulic` and drive loss `electrical - shaft` delivered to the
well-mixed finite reservoir. The complementary fractions leave to the explicitly
unmodelled surroundings. The operating point reports all input powers, both
losses and both heat destinations. The finite-reservoir implicit solve includes
this heat at every current material/flow trial, so heating changes coolant
enthalpy, viscosity and the next pump operating point. Its report retains
`drive_heat_j` separately from user-supplied external heat and branch return heat.
Do not include this same heat again in the external-heat argument.

The spatial stored-solid/reservoir driver includes drive heat in its joint
energy balance. Hydraulic output still belongs to the original network's
hydraulic-dissipation or explicit omission model and is never added twice.
All existing source, work, energy and rollback gates apply. This is an
endpoint-held quasi-steady equipment model with declared heat routing, not a
rotor/inverter transient, cavitation model or automatic map extrapolation.
Focused pump, reservoir and stored-cooling tests exercise input conservation,
shutoff consumption, interior-curve infeasibility, thermal feedback, jointly
stored energy and transactional refusal.

## Finite two-stream heat exchange with original caloric sources

`conjugate::caloric::exchanger::CaloricHeatExchanger` solves steady co-current
or counterflow exchange between two `ExchangerStream` owners with prescribed
positive mass rates and total UA. Each retains its original
`IsobaricHeatCapacityLaw`, supported temperature interval, source knots and
fixed composition. The fluids may have different isobars and independent
enthalpy datums: they exchange heat, not mass.

One shared heat rate determines both original relative-enthalpy outlet
inverses. The required conductance is the integral of heat progress divided
by the actual local stream-temperature difference. Counterflow reverses cold
progress along the hot direction. The integral is split at BOTH original
Cp-knot sets. On each paired affine-Cp interval, equality of the positive
capacity rates locates every possible interior stationary temperature gap;
positive inlet/outlet gaps alone cannot admit an interior thermal pinch.

Constant Cp uses the stable analytic NTU relation, including equal and nearly
equal capacity rates. Variable Cp uses a bounded scalar bracket and Simpson
quadrature with analytic fourth-derivative truncation bounds plus an explicit
floating-point allowance. Heat and conductance acceptance gates both apply.
The numerical allowance is not an interval-arithmetic certificate or source
uncertainty. Quadrature work and nonlinear iterations have explicit caps,
with cancellation during source preparation, panel work and publication.

The report retains the reciprocal signed heat, both relative enthalpy changes,
rounded outlet temperatures, original-source endpoint caloric residuals and
their explicit numerical allowances. Original entropy integrals include the
stated first-order correction for unresolved endpoint energy, with a separate
signed entropy allowance; negative production is never clipped. Endpoint
rounding allowances do not relax the independent Q/UA gates. Original source
support limits both streams, and further exchange outside either law refuses.

This is a held-UA, steady, single-phase exchanger. It introduces no wall/fluid
inventory, pressure loss, pump work, phase change, external heater or mass
mixing. Extreme constant-Cp saturation that rounds a terminal temperature gap
to zero remains an explicit numerical refusal because the represented heat
cannot support the required finite conductance inversion.


Nine focused regression tests cover independent constant-Cp NTU limits, equal
and nearly equal capacity rates, reciprocal stream reversal, weak heat below
absolute-temperature resolution, nonlinear manufactured heat and entropy,
actual original-Cp slope changes, distinct datums and isobars, interior
counterflow pinch, both source-domain limits, zero transfer and explicit
saturation, bounded work, nonconvergence and cancellation. The three Rust
sources were parsed with the pinned formatter. Native Cargo compilation and
test execution were not available in this resource-exhausted session.

## Finite-throughput original water/steam condensers and boilers

`conjugate::caloric::exchanger::phase::WaterPhaseHeatExchanger` couples one
original nonlinear-Cp `ExchangerStream` to a finite positive water/steam mass
flow with prescribed finite nonnegative UA. The water inlet must be an
original saturation state, including either pure saturated endpoint, from
the Gibbs, dense Helmholtz or regional IF97 owner. Its same original pressure
chart resolves the outlet through the unified p,h inverse; the saturation
temperature stays fixed while quality changes.

The sensible path reuses `IsobaricHeatCapacityLaw::linear_heat_exchange`
directly, with exposure `UA / sensible_mass_rate`. All original affine-Cp
segments and source boundaries retain that owner's native integration and
bounded inversion. No new caloric quadrature, EOS, latent Cp or temperature
table is introduced. Co-current and counterflow have the same solution while
the entire water path is isothermal. The two rates and the UA are finite;
water is not an inexhaustible reservoir. The retained water enthalpy request
is compared to both original saturated endpoints BEFORE rounding absolute h.
Requests that exhaust available liquid or vapor refuse transactionally.
There is no automatic clamp, subcooling zone or superheating zone.

One signed Q, positive from water to the sensible stream, supplies equal and
opposite relative enthalpy changes. The report separately audits original
represented endpoint h and water u+p*v, so pressure-flow work is included
once. Numerical endpoint rounding is exposed in watts. Each original
requested-pressure closure defect has a separate caller budget in J/kg;
their mass-weighted absolute sum is the stated additional allowance on the
independent water first-law check. The original pressure state is retained,
without an energy projection.

The signed raw water entropy change uses the original plateau delta_s/delta_h
and retained Q, preserving small changes before absolute endpoint rounding.
The sensible side reuses the original Cp/T entropy integral and its existing
endpoint correction. Raw total entropy production has an explicit caller
allowance in W/K, defaulting to zero beyond numerical roundoff. The separate
local thermal audit uses water delta_h/T_sat and reports the corresponding
original saturation seam adjustment. This pressure-specific primitive never
becomes a global source entropy or erases a negative raw balance. Neither
entropy report is clipped. Source inverse work counts and cancellation are
retained, and all gates precede result publication.

This is steady heat exchange with prescribed quality-independent UA; it does
not infer a boiling/condensation coefficient, phase slip, pressure loss,
wall storage or mass mixing. Focused regressions cover independent NTU and
original latent-energy predictions in both heat directions, nonlinear Cp
knot crossings with independent datums, dense and regional source retention,
latent exhaustion below absolute-h resolution, signed source-entropy seam
admission, finite source support, zero transfer and cancellation.

## Continuous co-current water/steam exchange across phase boundaries

`conjugate::caloric::exchanger::co_current_water::CoCurrentWaterHeatExchanger`
accepts an original Gibbs `WaterState`, its finite positive mass rate, an
original-Cp `ExchangerStream`, and prescribed total UA. Water may start and
finish in different sensible or latent branches. Desuperheating, condensation
and subcooling, or feedwater heating, boiling and superheating, all use the
same original enthalpy coordinate and source domain.

For normalized cumulative conductance x from zero to one, the steady
co-current equations are `m_dot_water dh_water/dx = -UA*(T_water-T_sensible)`
and the equal opposite sensible-stream equation. These are exactly the
existing conservative two-body thermal equations under the declared map
`m_aux=m_dot*tau`, `t_aux=x*tau`, with `tau=1 s`. The auxiliary time and mass
are numerical coordinates; they do not prescribe exchanger inventory,
physical transit time, axial mixing or stored energy. Native unit auxiliary
geometry supplies topology only. All external boundary heat is zero and the
only active link has the caller's total UA.

The original heterogeneous `DepositionSimulation` and native adaptive SSPRK2
owner execute every stage. Both original material `ThermalLaw` implementations
resolve their own T(h), with the actual water plateau and original Cp knots.
The pairwise step limit additionally uses the SUM of the two original global
inverse-enthalpy slopes, preserving passive co-current progress through every
Euler stage. Native coarse/fine enthalpy and optional temperature norms,
datum-independent relative accuracy, minimum progress and shared trial/stage
limits remain explicit. The accepted conservative fine trajectory is used;
no extrapolated state, separate quadrature or new energy integrator is added.
These are local numerical error estimates. Phase corners require refinement;
the report does not promote them into a certified global heat error.

Original water `advance_heat` independently audits stable source and represented
endpoint changes in h and u+p*v. The original sensible endpoint audit checks
the same reciprocal heat. Original raw water entropy, its local fixed-pressure
thermal seam, and the sensible Cp/T integral retain separate signed reports.
Endpoint inverse/rounding allowances are explicit, and the caller's raw-source
entropy allowance is not replaced by a seam correction. Every gate and final
cancellation checkpoint precedes result publication.

Source-domain exhaustion, unrepresentable nonzero heat, passive-gap crossing
and exhausted work refuse. This owner currently admits the continuous original
Gibbs chart; dense/regional extensions require their own continuous-path
admission. It models co-current plug flow with prescribed UA and fixed
composition/pressure on each side. Pressure loss, wall inventory
and boiling/condensation correlations are separate capabilities. The existing
plateau exchanger remains the exact native fast path for wholly isothermal
water exchange. Focused regressions compare that plateau solution, independently
integrated full-path UA duties, accuracy refinement, mass-rate/UA scaling,
independent Cp datums, source limits, zero exchange and cancellation.

## Continuous countercurrent water/steam exchange across phase boundaries

`conjugate::caloric::exchanger::countercurrent_water::CountercurrentWaterHeatExchanger`
accepts the same original Gibbs water inlet, finite water mass rate and
original-Cp sensible stream with independent fixed pressures and temperature
domains. `solve` finds the reciprocal duty for prescribed total UA;
`rate_duty` reports required UA for a specified signed water-to-sensible duty.
Both admit desuperheating/condensing/subcooling and feedwater/boiling/superheating
without clipping exhausted phase inventory or replacing latent heat by Cp.
Countercurrent outlet temperature crossover is physical and is admitted when
the entire local paired-stream gap remains positive.

With duty magnitude Q and sign sigma, the water-direction coordinate q runs
from zero to Q. Original water enthalpy is `h_w,in-sigma*q/mdot_w`; the opposing
stream has remaining native heat progress `Q-q`. Required conductance is
`integral(dq / (sigma*(T_w(q)-T_s(Q-q))))`. Native water p,h inversion and the
existing Cp `SourcePath` retain their original source support. Every original
saturation endpoint and Cp knot splits the integral. Both oriented source
temperatures decrease with q, so cross-endpoint temperature ranges on each
accepted quarter-panel bound its WHOLE local gap. Endpoint-only pinch tests
cannot admit an unsampled interior crossing. Unresolved ranges subdivide;
sampled pinches, inverse-resolution limits and exhausted work refuse.

The original exchanger's fixed-order adaptive traversal is shared, while its
existing affine-Cp panel/error formulas remain unchanged. Water panels use
embedded composite-Simpson error ESTIMATES plus a separately reported allowance
from native inverse residuals, original global inverse-temperature slopes and
rounding. These are not analytic truncation bounds or interval certificates.
The minimum admitted temperature gap and numerical duty bracket retain that
qualification. Monotone original-source ranges supply the early lower-UA
rejection; estimated partial quadrature does not. The monotone scalar duty
solve uses both UA and heat-width gates, original finite-throughput brackets,
per-UA panel/outer-iteration bounds and one water-query budget shared across
preparation and every accepted or rejected trial.

The existing original endpoint audit is shared as `WaterExchangerState`;
`CoCurrentWaterExchangerState` remains its public compatibility name. Original
water `advance_heat` checks stable and represented h and u+p*v changes, the
original Cp endpoint checks the reciprocal duty, and raw source entropy and
thermal seam correction remain separate signed accounts. Pressure flow work
is included once. Countercurrent results have no SSPRK2 trajectory. Final
native calorimetry retains its own bounded inverses and cancellation, and no
partial result escapes any failed gate.

Focused tests independently integrate full condensation and boiling duties,
exercise physical outlet crossover, compare the exact plateau owner, preserve
mass-rate/UA scaling and independent Cp datums, cross nonlinear Cp knots with
distinct source domains, reject a positive-endpoint interior pinch, and check
source exhaustion, query limits, zero exchange and cancellation. The same
Gibbs-chart, prescribed-UA, fixed-composition and fixed-pressure limits as the
co-current owner apply; no dense/regional seam continuation, pressure-loss law,
wall storage or empirical boiling/condensation correlation is inferred.

## Prescribed mass-flow convection on native conducting surfaces

`implicit::PrescribedMassFlowConvection` borrows an original prepared gas or
liquid exterior owner and accepts a nonnegative prescribed mass-rate magnitude
plus the stream's positive cross section. Every actual bath-temperature query
resolves the SAME source at that bulk temperature and uses
`U=mdot/(rho_bulk A_flow)`. Film conventions still use the original film density,
viscosity, conductivity and Cp in Re/Pr/h; film density never determines the
bulk stream velocity. The flow cross section is distinct from the native
heated area. Bulk-held laws reuse the already queried snapshot.

Existing forced pointwise, sphere, plate-strip and laminar/turbulent circular
section laws retain their exact geometry, entrances, source support, thermal
direction, wall-viscosity correction and empirical domain. Negative rates and
natural/mixed geometry refuse; zero flow gets no speed or convection floor.
The magnitude does not infer a reversed leading edge. The bath, mass rate and
flow area remain fixed in the native surface derivative, so dU/dTs=0 and the
existing source-owned film/wall chain is retained. No bath or hydraulic-network
Jacobian, pressure solve, fan power, material replacement, finite fluid store or
new heat integrator is introduced.

The diagnostic response includes the actual bulk snapshot, speed, original
coefficient evidence/surface tangent and optional analytic dh/dmdot at fixed
temperatures and geometry. An unavailable original Re derivative remains
`None`, including source-domain endpoints, without blocking an independently
valid held-rate surface solve. Finite sections use their own coefficient
partial; Hausen uses the original correlation owner's stable section tangent.
The native conducting surface continues to own solid/coating resistance,
radiation, caloric identity, cancellation and surface balance.

Four focused regression definitions cover bulk/film separation, actual
bulk-temperature updates, ideal-gas pressure cancellation at fixed mass flux,
gas/liquid surface tangents, all forced geometry rate partials, source/geometry
and zero-flow refusal, unavailable rate derivatives at a fixed regime endpoint,
and an original conducting two-bath surface against an independent integrated
conductivity/plate balance. The new Rust files were formatted and the touched
modules parsed; native compilation and execution remain blocked by the exhausted
8 GiB build environment. No executed test pass is claimed.

## Analytic two-stream inlet and conductance response

`CaloricHeatExchanger::solve_with_inlet_response` first admits the original
steady heat exchange, including source support, both Q/UA gates, endpoint
calorimetry and entropy. It then returns `ExchangerInletResponse` for signed
heat versus each inlet's specific enthalpy and versus held UA. This response
keeps the independent original caloric laws and datums; it does not differentiate
an averaged Cp or an outlet-temperature chord.

For varying Cp, bounded quadrature integrates inverse-capacity temperature-gap
moments on the same original source breaks and interior-pinch geometry as the
primal solve. Differentiating the implicit UA equation includes the cold-side
counterflow dependence on total heat. Stable integral identities retain both
own-stream outlet-enthalpy gains directly, including a weak finite-conductance
anchor that would disappear in subtraction from a rounded unit transfer gain.
Constant Cp, zero UA and equal inlet temperatures have their analytic limits;
equal-temperature response uses the actual local Cp of each source.

`ExchangerResponseConfig` controls the absolute and relative enthalpy-gain
accuracy. Analytic derivative bounds and explicit floating-point allowances
propagate into the reported heat sensitivities, both own-stream gains and
conductance sensitivity. Cancellation and original quadrature/work bounds
remain active. Insufficient accuracy or unresolved arithmetic refuses the
response. The derivatives are evaluated at the admitted heat: the separately
reported primal heat error is not silently incorporated as a certified
sensitivity interval. These allowances are numerical estimates with analytic
truncation bounds, not interval arithmetic or material-source uncertainty.

The ordinary `solve` path shares only the source-break and panel-geometry
helpers and does not compute derivative moments. Six added regressions cover
independent NTU derivatives and reciprocal balances, independently integrated
nonlinear moments, stream reversal and independent datums, actual Cp kinks,
equal-temperature and near-equal-capacity limits, weak transfer and retained
counterflow anchors, zero UA, accuracy refusal and cancellation. All four
changed Rust files were parsed with the pinned formatter. Native Cargo
compilation and test execution were unavailable in this resource-exhausted
session.

## Finite exhaust heat recovery on native source-caloric networks

`CaloricTransportNetwork::with_heat_recovery` binds
`CaloricHeatRecoveryLink` values between whole external exhaust and fresh
supply ports of the existing directed-acyclic hydraulic graph. Actual external
MASS flows feed the original two-stream exchanger. Each port has at most one
owner, ordering is deterministic, and the native 64-supply reduced-solve cap
applies. Finite heat recovery, imposed mass returns and closed circulation are
exclusive modes in either binding order.

Each trial retains the network's one original isobaric Cp law, relative
enthalpy coordinates, native branch heat exchange and mass-weighted enthalpy
mixing. One signed exchanger Q debits the exhaust and credits its paired fresh
supply, including reverse transfer when the fresh stream is hotter. The
conditioned supply enthalpies close a reduced damped Newton system using the
analytic exchanger inlet response and the original DAG branch gains. The
shared inlet-response traversal preserves the existing recirculation arithmetic.
Retained own-stream exchanger gains and double-double matrix assembly keep a
small finite-UA feedback anchor from disappearing in subtraction.

`CaloricHeatRecoveryConfig` bounds accepted iterates, damped trial counts and
original exchanger work. A wall-path seed refusal may use one bounded native
zero-wall-heat recovery solve to obtain a supported initial supply. Final
acceptance requires both the sum of absolute supply-closure watts and a local
Newton-correction temperature test using the positive original minimum Cp.
The latter is a local convergence measure, not a global state-error proof.
Source support, exchanger pinch and Q/UA accuracy, entropy, cancellation and
original native branch admission all remain active.

`TransportMarch::heat_recovery` supplies a distinct
`CaloricHeatRecoveryReport`. The unchanged march fields retain the INNER
network balance with conditioned supplies and pre-exchanger exhausts.
The report independently recomputes the OUTER balance from original fresh
inlets and retained post-transfer exhaust enthalpies. Both physical watt
gates must pass; no exchanger numerical allowance enlarges either budget.
The report retains each exchanger's original calorimetry, response, entropy
and numerical diagnostics, actual supply closure, and a separate exhaust
temperature/enthalpy handoff rounding audit. Summed exchanger entropy is
explicitly limited to the recuperators; it is not a complete wall/junction
entropy account.

The existing `march`, initialization and shared-solid `solve` paths use this
composition, including the native coupling/acceleration driver. The graph's
hydraulic flow, density and branch film coefficients remain prescribed. This
held-UA steady recovery owner adds no mass-return edge, pressure-loss model,
fan work, changing composition, phase change, reservoir or exchanger wall/fluid
storage. Focused regressions cover independent finite-NTU closure, nonlinear
original Cp and datum changes, unequal mass rates, multiple port ownership,
shared-solid use, reverse transfer and transactional refusal. Changed Rust
sources were parsed with the pinned formatter; native Cargo compilation and
test execution were unavailable in this resource-exhausted session.

## Pressure-selected heated gas with source density

`medium::steady_gas::pressure::solve_pressure_driven_material_gas` selects a
positive mass rate on a constant-area axial path without side injections.
Every flow trial reuses the original stationary species/NASA9 enthalpy,
convection and resistance solve. Its retained local EOS density determines
velocity and the momentum term `beta * (mdot/A) * (u_out - u_in)`; this joins
the actual hydraulic loss sum in the imposed static-pressure balance.
The returned native gas field and transport observations are the accepted trial,
with no extra solve or material reconstruction after pressure convergence.

The caller authors one constant velocity-profile shape through its momentum
factor beta and kinetic-energy factor alpha. Both must be finite and at least
one, with `alpha >= beta²` up to 64-epsilon relative roundoff for unidirectional
profiles. Uniform flow uses (1, 1); a parabolic circular profile uses (4/3, 2).
Heat-transfer correlation names never select these factors. The declared factors
are retained in the result and apply to momentum and kinetic power consistently.

The existing `PressureDrivenBudget` supplies the positive bracket, pressure
tolerances, solve limit and hydraulic-omission allowance. A separate total
native-substep cap covers all flow and material-feedback trials. The source's
common thermodynamic pressure and native low-Mach bounds remain active;
`sum(dp_loss) + max |beta*G*(u_cell-u_in)|` limits pressure excursions without
cancellation. The sum of absolute regional kinetic-power changes has its own
omission allowance, so heating/cooling cancellation cannot hide omitted energy.
Every evaluated trial must be admitted; failures and cancellation preserve inputs.

Bisection selects an admitted sign-changing bracket without claiming monotonicity,
uniqueness or stability. This is a reduced common-pressure operating point, not
a compressible pressure profile, choking model, resolved velocity/thermal boundary
layer or complete energy solve. Existing tests add independent Poiseuille density
scaling, heated/cooled mass/enthalpy/momentum balances, retained profile factors,
and bracket, profile, shared-work, pressure, power and cancellation refusals.

## Reynolds-dependent turbulent resistance in native pressure graphs

`LossGraph::new_mixed` accepts `PassiveBranch::Turbulent`, retaining a
`TurbulentPipeLoss` with the original `MaterialFlowResistance::circular_turbulent`
geometry, equivalent sand roughness, Reynolds support, coherent material
snapshot and provenance. Its signed Darcy-Weisbach pressure loss and passive
`dp Q` use the same original Colebrook equation and fixed-material tangent.
Source-bound edges share one conserved momentum density; gravity must use that
same reference. The material-response compiler retains this nonlinear branch
when viscosity changes, without fitting constant linear/quadratic coefficients.

The existing nodal bracket solver remains the network owner. Algebraic
inversion of the original Colebrook equation avoids an inner friction root
solve at each pressure trial. Private monotone continuations below and above
the authored Reynolds interval support numerical bracketing through zero;
they are not laminar, transition or extrapolated physical models. Every
accepted graph/fan solution re-evaluates the original Reynolds support,
pressure, positive differential resistance and dissipated power. Unsupported
roots refuse, including zero/transition flow. Exact endpoints are rounded
inward and inverse roundoff is bounded before the unchanged source admission.

Fan intersections use actual nonlinear graph pressure trials rather than
quadratic homogeneity. Native tangent/adjoint solves reuse the existing
grounded Laplacian with `dQ/d(dp) = 1/(rho_ref d(dp)/dm)`. A resistance-scale
direction multiplies the whole pipe loss, such as changing length at fixed
diameter, roughness and fluid; it is not a viscosity or roughness derivative.
The thermal circulation closure compares fresh original pressure and
differential resistance at actual solved flow, in addition to its unchanged
coefficient and hA gates. Its existing omitted-dissipation watt budget still
accounts for actual `dp Q`; no power is silently deposited twice.

Constant-coefficient axial aggregation still refuses nonlinear pipe branches.
The original inertial graph also accepts these pipes through its native
midpoint/secant closure, as described below. This feature adds no
second network solver, transient turbulence model, rough-wall heat-transfer
correlation, source-domain extension or uncertainty certificate. Five added
regressions in existing test files cover junction/reversed flow and pressure
work, original length/pressure derivatives and adjoints, viscosity refresh,
fan closure, original Reynolds endpoints/refusals, coefficient-only owner
refusal and nonlinear thermal-response confirmation. The touched Rust sources
were parsed with the pinned formatter; native Cargo test execution was not
available during this batch.

## Source-density stack and circulation head on native gas flow

`solve_pressure_gravity_driven_material_gas` uses the same pressure iteration,
conservative thermochemical transport, source convection and retained velocity
profile factors as `solve_pressure_driven_material_gas`. The original entry
still requires positive imposed pressure; gravity explicitly admits signed or
zero imposed pressure with the same positive-flow bracket.

One authored `SeriesGravity` segment maps to each native axial cell. Segment
endpoint displacement cannot exceed the cell path length. The original gravity
owner evaluates `sum(rho_cell * g dot dr)` using accepted local bulk source
densities, with no reference-density buoyancy approximation. Pressure closure is
loss plus thermal acceleration equals imposed pressure plus hydrostatic gain.
A chimney's outside-column pressure difference is supplied explicitly; no ambient
density, vertical axis, gravity magnitude or inlet/outlet pressure is invented.

The common-pressure envelope adds absolute regional hydrostatic gains, preventing
nearly cancelling stack or loop heads from hiding large local pressure changes.
The retained original hydrostatic result reports conservative endpoint gravity
work. Potential-energy flow is its negative; an independent explicit watt
allowance bounds the sum of absolute regional potential-flow changes. Gravity
work is never added as heat. A closed geometric path has exactly zero net gravity
work but can have source-density buoyant head and nonzero regional omission.
The new read-only `section_specific_work_j_kg` accessor exposes the gravity
owner's original prepared regional work without repeating geometry arithmetic.

All native low-Mach, common-pressure, caloric/source, Reynolds and work limits
remain active. Every pressure trial must meet the gravity/energy allowances;
refusals and cancellation leave inputs unchanged. This remains a prescribed-inlet
reduced steady path, including on closed geometry; it does not enforce periodic
thermal state, evolve a finite inventory, resolve a continuous hydrostatic EOS
profile or close the omitted mechanical terms in total enthalpy. Refine the
native cells for density gradients. The existing gas battery adds zero-gravity
parity, independent downhill Poiseuille/density/signed-drive oracles, a sourced
heated stack and neutral-column refusal, plus closed-loop buoyancy, regional
energy/pressure admission, geometry and cancellation checks.

## Parallel source-density gas branches at common pressure ports

`steady_gas::pressure::parallel::solve_parallel_pressure_driven_material_gas`
selects positive mass flow independently in each native axial branch at one
common imposed static-pressure difference. The branch pressure/gravity owners
retain local source density, viscosity, convection, thermal acceleration and
the original enthalpy/species solution. Different temperatures or areas change
the actual split; no frozen volume-flow fraction is supplied. Zero branches
refuse and a single branch is supported. The complete set is bounded by
`MAX_PARALLEL_GAS_BRANCHES` and the native total steady-cell cap.

Branches must retain identical source objects, external conditions, common EOS
pressure, phase and inlet composition/temperature. Original branch pressure
controls are admitted before shared-work clipping. All branches either omit
gravity or author paths with identical acceleration and physical inlet/outlet
vertices; interior routes can differ within native cell-length admission.
Signed or zero drive remains confined to the original gravity entry. Every
branch receives only the remaining aggregate full-solve/substep allowance.
Any failure terminates the calculation without partial publication or retries;
accepted native solutions are retained without an additional flow solve.

The ideal adiabatic, zero-storage outlet mixes the same source's specific
enthalpy. `GasEnthalpyLaw` integrates original NASA Cp for each inlet-to-outlet
increment before mass weighting. The existing pressure/enthalpy-increment
reconstructor finds the common mixture and retains gas/water phase admission.
Formation-energy origins are not subtracted after rounding, and neither
temperature nor volume flow substitutes for conserved mass/enthalpy.
Reported inlet and mixed-outlet volume rates use their actual source densities.

Aggregate omission limits sum hydraulic dissipation, absolute regional kinetic
changes, and absolute regional potential changes over all accepted branches.
The kinetic bound additionally includes each branch's inlet and outlet kinetic
energy flux magnitudes, using its original profile factor. This explicitly
bounds energy ignored by the ideal split/merge; it is not a hydraulic junction
loss, manifold momentum closure, or kinetic heat deposition. Any authored
hydraulic junction loss remains in the native branch resistance. These are
common-pressure, low-Mach positive-flow ports, not a general pressure graph,
recirculating inventory, compressible pressure field, or complete manifold
momentum/energy solve.

The existing gas battery adds independent parallel Poiseuille and density
scaling, scalar parity, nonlinear NASA enthalpy mixing, thermally redistributed
flow, common-port gravity, discarded manifold kinetic flux, shared-work and
source/geometry/control/cancellation refusals. Product and tests were parsed
with the pinned formatter; native Cargo compilation/test execution was not
available in this resource-exhausted session.

### Source-bound fan operating points on heated parallel gas branches

The `steady_gas::pressure::parallel::fan` owner now connects the existing native
fan bank to the original parallel gas pressure/thermal solve. Each trial imposes
one positive static-pressure difference on the same branch ports, resolves the
original thermochemical fields and flow split, and binds the reference bank with
`FrozenFlowMedium::bind_fan_bank` to the **exact resolved network inlet**. Actual
inlet volume flow selects the fan point; neither mixed-outlet volume nor mass
divided by reference density is substituted. Original fan curve provenance,
arrangement, speed, stall support, pressure tolerance, density/Reynolds similarity
and pressure-fraction admission survive. The accepted native network is returned
directly, retaining its source enthalpy mix, gravity and omitted-mechanical-work
accounting. This remains the common-pressure, low-Mach, constant-area-per-branch
model with ideal imposed-static-pressure manifolds, not a manifold momentum or
hydraulic junction-loss closure.

The caller supplies a finite positive pressure sign bracket, explicit residual
tolerances and bounded complete network evaluations. Native flow evaluations and
thermal substeps share one outer allowance across every pressure trial and
branch; original aggregate controls are admitted before remaining work is clipped.
An unsupported or failed trial stops without extrapolation, partial publication
or a second accepted-state solve. Heating and buoyancy can make the response
nonmonotone: the result certifies neither uniqueness/stability nor propagation of
the retained fan-curve uncertainty.

A static pressure curve does not provide shaft efficiency or a fan heat model.
The caller must author an independently applicable minimum **bank static
efficiency** in (0,1], valid at every evaluated flow/speed/medium, and a maximum
neglected shaft-power allowance. At actual inlet volume Q the declared upper
power is Q times nominal static rise times (1 + retained relative pressure
tolerance), divided by this minimum efficiency. This conditional allowance is
reported separately from branch/manifold hydraulic, kinetic and potential
omission bounds; it is never deposited as heat or subtracted again. The source
inlet thermodynamic state stays frozen across the admitted omitted fan heating
and compression. Zero static supply, an absent/invalid efficiency bound or an
exceeded shaft-power allowance refuses; motor electrical losses are not inferred.

Focused regressions cover an independent fan/Poiseuille intersection under
density and speed changes, heated branch redistribution with nonlinear original
enthalpy, actual inlet-volume selection, native source/stall refusals, and shared
work, cancellation and separate omitted-fan-power admission.

## Explicit continuous circular-pipe viscosity law

`MaterialFlowResistance::circular_churchill` selects Churchill's 1977 empirical
laminar/transition/turbulent correlation with authored Reynolds support inside
[0, 1e8] and equivalent sand roughness / diameter inside [0, 0.05]. Existing
Poiseuille and Colebrook constructors retain their original applicability.
This assumes steady developed Newtonian pipe flow; it does not resolve
transition intermittency, entrance profiles, compressibility or unsteady stress.
The transition fit is an engineering approximation, not an onset criterion.

The implementation scales the exact Poiseuille geometry by `f_D Re / 64`,
which has a finite unit limit at rest. Original source dynamic viscosity and
momentum density determine Reynolds number, signed pressure and pressure work.
The exact fixed-material derivative differentiates the same correlation,
including its Reynolds dependence. Weak flow and rest recover the original
Poiseuille resistance without artificial drag or a singular friction factor.

Existing series and native graph owners retain this original nonlinear law and
its material snapshot. Series inversion intersects authored Reynolds intervals;
graph inversion uses bounded Newton/bisection inside admitted flow endpoints.
Reverse flow remains dissipative. An explicitly admitted rest state now passes
native pressure/work closure and graph sensitivity uses its positive original
differential resistance. Constant-coefficient exports continue to refuse
nonlinear friction; a hydraulic regime selection does not extend any separate
heat-transfer correlation's domain.

Regressions compare the independent published Darcy equation, perturb the
analytical tangent through transition, check rest and signed-flow symmetry,
and exercise actual viscosity changes, source-retaining series pressure solves,
native graph inversion and the zero-flow Poiseuille sensitivity.

Reference: S. W. Churchill, "Friction-factor equation spans all fluid-flow
regimes" (1977), equations (5), (18), converted from the paper's friction
convention to Darcy friction:
<https://files.engineering.com/files/85c0f3a6-a102-4a22-9d35-f15858c0dd2b/CEM_-_Friction-factor_equation_%281977%29.pdf>.

## Liquid wall properties in native developing graph convection

`DevelopingDuct::TurbulentCircularLiquidWall` explicitly selects Gnielinski's `(Pr_bulk/Pr_wall)^0.11` correction for the original turbulent circular section. The original `TurbulentCircular` path retains K=1. Only the native liquid graph owner admits the new selection; gas names, compositions and numerical Prandtl values never establish liquid eligibility.

Both prescribed-wall and finite-wall paths keep the original bulk material state for Re, Pr, conductivity, density admission and the accepted native flow. The same fixed-pressure/composition liquid owner is additionally queried at the actual fluid-facing wall temperature. Its retained source snapshot supplies wall Prandtl number, including Cp and conductivity, rather than substituting wall viscosity alone. This query introduces no mass, caloric store, hydraulic density replacement or independent heat source.

Every original signed flow endpoint still admits the same physical entrance distance and corrected section; a trajectory crossing zero cannot hide outside the forced-flow domain. Full wetted-area geometry and the original source and correlation bounds remain active. Samples retain both cumulative K=1 means, the actual corrected section and coefficient evidence, the actual wall source, and the correction's independent analytic input partials. Those partials are not total derivatives through the material or graph.

The existing conjugate normal-path solve re-queries this same corrected section at each actual surface trial and publishes only its combined reciprocal conductance to the native wall/fluid link. The bare film hA remains diagnostic. Locally changing bulk or wall properties remains a held-section approximation; it does not resolve axial wall spreading or a nonisothermal boundary-layer history, and corrected model discrepancy remains unbounded.

Focused regressions cover independently evaluated section heat, original wall-source coordinates, signed entrances, preserved K=1 behavior, source support and gas refusal, accepted native momentum/advection/storage balance, and conjugate surface/reciprocal heat balance in heating, cooling and equilibrium. They do not establish physical calibration.

## Prescribed stream mass flow on finite gas convection

`FilmGasSurface::with_prescribed_mass_flow` binds a positive stream or recirculation mass rate and flow cross section to an existing forced heat law. Every nonlinear trial resolves velocity as `mdot/(rho_bulk A_flow)` from the same sealed original gas state's actual EOS density. It does not query unused bulk transport merely to infer density. External-film correlations retain their separate film density, viscosity, Cp and conductivity. Internal circular and rectangular flow areas must agree with the supplied passage geometry within roundoff; physical heated entrance coordinates remain unchanged. Natural and direction-dependent mixed laws refuse this unsigned forced-speed selection.

The returned sample retains the prescribed inputs, actual bulk density and velocity. The default authored-speed API and shared piston/body and fixed-volume wall-feedback algorithms remain intact, with the original source checks, regime bounds and exact zero thermal-area insulation. Prescribed rate controls convection only: it supplies no open mass or enthalpy flux, fan/pump work, axial solution or new energy owner. Focused regressions cover independent film transfer, mass-flux invariance, passage geometry, original variable-volume work/energy, and fixed-volume energy/entropy conservation. This does not close the broader conjugate-channel validation bead.


## Source-resolved dilute forced spherical evaporation

`mass_transfer::sphere::ForcedSphereMassTransfer::Frossling1938` explicitly
selects the original direct mass-transfer wind correlation. Its diameter-based
form is `Sh = 2 + 0.552 sqrt(Re) Sc^(1/3)`, with `k_m = Sh D_AB/d`.
This rung admits Re in [2,800] and Sc in [0.6,2.6]; it never switches models,
clamps unsupported states or substitutes a heat-transfer analogy. The scalar
correlation performs bounded work without allocation. Area-averaged Sherwood
and finite kinetics form a lumped series closure; angular variation is unresolved. The original paper and
chosen range are documented beside the implementation.

`evaporation::droplet::forced::SourceForcedSphericalDropletFilm` borrows the
same pinned binary diffusivity source as the stagnant sphere. It resolves
molecular D and original gas density/viscosity at the arithmetic mean of actual
surface and bulk temperature, with actual bulk pressure/composition. The mean
is also the representative Stefan-film temperature; actual interface
temperature, finite accommodation kinetics and both caloric owners stay
original. True D, explicit Sh, transport groups and original material receipts
remain separate. No parcel is heated merely to evaluate its film properties.

Relative speed and a positive dilute-vapor bound are explicitly authored.
The bound cannot exceed 0.05 and applies to bulk vapor and equilibrium surface
pressure fraction at every phase-film evaluation. This ceiling is an
implementation policy, not an experimental accuracy certificate. Initial and
accepted-endpoint source/correlation support are required; a late refusal
discards the complete original phase/heat candidate. Source claims never
silently switch. Optional sensible heat remains independently authored.

Initial geometry and transport freeze for one first-order interval; final
samples describe the next interval. This is an isolated, uniform, quasi-steady
sphere in dilute binary continuum gas with nearly constant film properties
and negligible buoyancy, deformation and compressibility. It supplies no
trajectory, drag, strong-blowing correction, breakup, nucleation or gas-volume
work. Focused numerical/source regressions cover transport scaling and range,
mean-film sources, both phase directions with finite kinetics, reciprocal
heat, late refusals and original IF97/NASA calorics. The water fixture now
constructs mole fractions by original species identity. All new regression
execution remains unclaimed because the execution workspace is unavailable.

## Finite-cell liquid wall-Pr convection

Source-backed deposition boundaries and reciprocal finite-fluid links can select
`SourceCellConvectionLaw::ForcedTubeSectionLiquidWall` explicitly. The existing
Gnielinski turbulent circular-section owner receives bulk Re, Pr and conductivity
from the original fluid source, while that same resolver supplies wall Pr at the
actual trial wall temperature. The correction is `(Pr_bulk/Pr_wall)^0.11`; the
physical inlet, original cumulative means, corrected section Nu/coefficient and
correlation evidence remain available. Unselected `ForcedTubeSection` retains
K=1 and its single bulk-property query.

Liquid eligibility comes from `ThermalFluidTransport::resolved_liquid_source`,
which defaults to no source and returns the retained `ResolvedLiquidState` for
the original liquid owner, including borrowed views. Admission verifies both
thermal snapshots, one explicitly typed absolute-temperature coordinate, one
immutable card, identical property requirements and viscosity basis, selected
claim receipts, and unchanged non-temperature coordinates. Numeric-only and gas
callbacks, inconsistent witnesses and switched source bundles refuse.

Boundary exchange supplies one external power; a finite link carries one
conductance with equal and opposite heat on the two existing caloric cells.
No surface fluid mass, second temperature history or solver is added. Prepared
implicit liquid boundaries retain their existing
`with_gnielinski_wall_correction()` selection; the callback-only enum variant
refuses at implicit construction. Original empirical and property-ratio admission
remains in force. Corrected section model discrepancy is unbounded; this adds no
single-phase stability claim or independent experimental accuracy bound.

Four focused regression definitions cover independent heating/cooling transfer,
original bulk/wall sources, both native finite stores and link reversal, gas and
numeric-source refusal, and changed source claims or non-temperature context.
Their execution is not claimed by this implementation checkpoint.

## Shared atmosphere for heterogeneous sourced spherical droplets

`evaporation::droplet::chamber` explicitly assigns each physical droplet its
existing stagnant or forced-sphere film adapter, original liquid law and
caloric inverse. Different sizes, temperatures, kinetics and prescribed
relative speeds share one original material chamber. Canonical stable IDs and
a caller-owned positive roster ceiling bound the work. Volatile identity is
admitted before the liquid source is queried.

Every initial film sees the same actual gas. Stagnant transport retains the
original bulk-temperature D/r path; forced transport reuses the pinned source,
arithmetic mean film temperature, original gas density/viscosity, selected
Sherwood correlation and existing dilute guard. That guard remains active at
every phase evaluation, including after collective sensible heat. One joint
material phase/heat owner handles shared vapor and caloric storage. Final
source/correlation support for every accepted size and common gas must pass
before any candidate is returned; no model or source claim silently switches.

Full typed initial/final film receipts retain molecular D and Sh separately.
Geometry and transport freeze through this first-order interval; endpoint
samples prepare the next interval. Gas volume is external and fixed. This
well-mixed humidity/heat feedback introduces no statistical weights, wakes,
shielding, trajectories, drag, coalescence or moving-volume work.

Five focused G0/G3/G4 regression groups specify mixed models/directions,
source-based rates, original scalar parity, shared reciprocal heat, canonical
ordering, early identity refusal and atomic late source/dilute refusal.
New Rust modules passed pinned rustfmt parsing/format checks; runtime and
full-crate compilation remain unverified in this environment.

## Prescribed mass flow with original tube and sphere wall properties

`ForcedTubeGasSurface`, its explicit `laminar()` selection, and
`ForcedSphereGasSurface` accept `with_prescribed_mass_flow` with a positive
stream mass rate and flow cross section. They share the film adapter's original
EOS-density primitive: every nonlinear trial resolves
`u = mdot/(rho_bulk A_flow)` before the same source is queried at the actual
body/wall temperature for viscosity. Bulk Re, Pr and conductivity remain bulk
properties. Tube flow area must match `pi D^2/4` within roundoff; a sphere's
external stream area is independently authored. Laminar whole-length heated
area admission remains unchanged.

Both selections reuse the original piston/body and fixed-volume wall-feedback
solvers, with original source, correlation and geometry admission, exact zero
thermal-area insulation, and gas mass/energy ownership. Samples retain rate,
area, actual bulk density and speed. Fixed-speed calls and the existing
`film::GasMassFlowState` path remain available. Prescribed circulation supplies
no open mass/enthalpy flux, fan work, axial flow solution, second heat integrator
or additional empirical gas validity for Sieder--Tate.

Four focused regression definitions compare independent Sieder--Tate/Whitaker
values and actual bulk/wall source coordinates, distinguish EOS density from
wall-property density, and cover both native storage owners, latent heat,
mechanical work, entropy and unchanged source/correlation refusals. Their
execution is not claimed by this implementation checkpoint.

## Original gas transport feeding a collective droplet atmosphere

`evaporation::droplet::chamber::transport` composes the existing original gas
network with the sourced collective droplet chamber in one immutable
candidate. The caller explicitly identifies one finite gas region. Every
authored transport path samples the original network together; all droplet
films then sample that region's actual transported temperature, pressure and
composition before one joint phase or heat-then-phase interval.

The accepted chamber gas moves into the complete transported parcel map,
retaining unrelated nodes and every original fixed volume. The original
transport owner rechecks its mechanism-specific final driving-force ordering
after droplet heat and phase change. No second flow is sampled or spent.
Late source, pressure-order, conservation or cancellation refusal publishes
no gas or liquid endpoint. Callers retry from the original states.

Original species enthalpy, latent and sensible transfers are applied exactly
once by their existing owners. The report combines signed residuals and
noncancelling local mass/species/energy defects under one caller-owned
whole-interval allowance; opposing stage errors cannot erase that budget.
Original substep receipts remain available. Transport's entropy receipt
describes its own substep, not the subsequent droplet process.

This is first-order transport then collective heat/phase splitting with
explicit interval refinement, source selection, relative speeds and heat
coefficients. It adds no trajectory, entrainment, drag, shielding, statistical
weights or variable-volume work. Three focused G0/G3/G4 regression groups
specify actual nozzle-fed source refresh and independent full inventories,
empty-path parity, and atomic late pressure-reversal refusal with retry.
New modules passed pinned rustfmt parsing/format checks; compilation and
runtime execution remain unverified in this environment.

### Source-conjugate finite-cell convection

`PreparedSourceConjugateLinks` binds original zero-conductance contacts to
source-qualified wall conductivity layers and their exact native caloric,
contact, area and geometry identities. The existing material integral-coordinate
solver eliminates the actual fluid-facing temperature against the existing
source-convection callback. Bulk properties remain at the finite-fluid
temperature; selected film properties use the surface/bulk mean, and wall
viscosity or explicit liquid wall Pr uses the solved surface. The repeatable
fluid callback must retain its native material and pressure/composition context.

Only the combined material/film conductance enters the native reciprocal link.
Retained diagnostics include the original conductivity line, solved surface,
fluid sources, surface-balance residual and actual native GΔT representation.
The original native caloric temperature remains authoritative even where layer
enthalpy admission snaps a source endpoint; an unsnapped trace outside original
conductivity support refuses without clamping. The existing graph fixed-film,
surface-feedback and developing-layer paths use the same native-trace rule.

`SourceConjugateHeatAssembly` supplies these links through the existing complete
`HeatAssembly`. Original finite stores, other heat mechanisms and coupled
drivers retain ownership. Existing prescribed-reservoir composition remains
available; finite G is not duplicated as external power or a source-rate bound.
Disjoint mechanical contact inputs pass through; a metric touching a conjugate
endpoint refuses until the corresponding physical geometry is supplied.

The wall remains homogeneous, source-free, pure-phase and quasi-steady. No new
solver, correlation, surface storage, contact resistance or lateral spreading
is inferred. Original source and correlation admission remains active.

### Conducting surfaces on the original finite gas/body owners

`FilmGasSurface::with_conducting_layer` and its prescribed-mass-flow counterpart
bind the original finite body's exact caloric object to a source-qualified
pure-phase conductivity layer and explicit planar or radial geometry. The
existing material kernel eliminates the actual fluid-facing temperature against
the original film. Only the combined conductance reaches the existing piston
or fixed-volume heat owner. Film properties, natural or mixed temperature
contrast, physical entrance coordinates, and bulk-density mass-flux conversion
retain their original owners.

Samples preserve the actual body and surface states, original conductivity
source, material/film power residual, and conductance representation error.
Native use checks the layer residual against the enclosing source budget in
joules. Initial, midpoint, and accepted endpoint source support remain
mandatory, including conductivity support at each actual native body
temperature when the layer's enthalpy admission snaps a source endpoint.
Zero film area retains prepared geometry/caloric binding and participating
body-state admission but bypasses the material root and gas-property/correlation
queries.

Isochoric use requires full positive exposed-area coverage or explicit
insulation, with body-node radiation disabled. The original Biot check uses
combined conductance per unit area, and existing entropy accounting describes
the complete passive material/film connection. This remains a quasi-steady,
source-free normal wall path; it introduces no extra thermal storage, radiation
relocation, wall spreading, transient wall field, independent contact-resistance
law, or additional heat integrator.

## Adaptive source-backed collective droplet evaporation

`evaporation::droplet::chamber::adaptive` connects the original stagnant and
forced spherical film adapters to the material chamber's common adaptive time
controller. Both phase-only and reciprocal sensible-heat entry points evolve a
fixed positive-droplet roster in one original fixed-volume gas parcel.

Each original trial interval receives its actual incoming droplet inventories
and gas state. Geometry, selected diffusivity, and the forced adapter's admitted
mean-film density, viscosity and Sherwood number are refreshed through the same
original source resolution. Source identities, finite interface kinetics,
correlation applicability and the actual-surface dilute guard remain mandatory;
no fitted replacement source or heat/mass-transfer analogy is introduced.

Only original positive phase step-size limits request smaller physical steps.
Source, domain and cancellation failures remain terminal with their original
error payload, including failures on a trial endpoint. The material controller
owns local time error and the complete-window physical account; rejected and
coarse trials contribute neither accepted transfer nor retained film receipts.
The final payload contains only the last accepted fine interval's original film
samples. A zero-duration request still admits original source and species
identity and returns its original samples without physical transfer.

Three regression groups cover exact two-fine-interval composition and refreshed
films, heated accuracy rejection/permutation/atomic cancellation retry, and
zero-time plus late source-domain refusal. Formatting and bounded static review
passed; compilation and runtime execution remain unverified because this
workspace has no runnable Cargo toolchain. The original source/correlation
limits and the material controller's local-error and positive-seed limits apply.

### Implicit source-conjugate finite-cell convection

`SourceConjugateHeatAssembly` enters the shared thermal-only implicit candidate
directly. Original wall k(T), the eliminated fluid-facing surface temperature,
source-pinned bulk/wall fluid properties and combined reciprocal G are
re-evaluated inside every native nonlinear trial and final acceptance.
The original finite wall/fluid caloric stores supply the energy unknowns;
bare-film hA remains diagnostic. Correlation validity, quasi-steady layer
assumptions and temporal/source refinement obligations remain unchanged.


### Specialized source films on conducting gas surfaces

Source-pinned finite-gas conducting surfaces accept the original turbulent
Sieder–Tate tube, explicit laminar whole-tube Sieder–Tate mean, and Whitaker
sphere selections through the sealed `ConjugateGasFilm` interface. Their
fixed-speed and prescribed-mass-flow forms reuse one `ConjugateGasSurface`
material elimination and the original piston/body or isochoric heat
transaction. Original bulk properties and source identities remain
authoritative; only the original wall-viscosity query moves to the eliminated
fluid-facing temperature. Whole heated length, laminar full-area admission,
sphere viscosity-ratio limits, and actual bulk-density velocity remain
unchanged. The existing `ConjugateFilmGasSurface` type and film sample payload
remain compatible for both legacy constructors. All shared original-caloric,
native-temperature conductivity support, pure-phase, radiation, Biot,
cancellation, insulation, and source-work admissions continue to apply.

### Thermal feedback into resolved pressure and head drives

`carreau::transient::adaptive::advance_resolved_with_thermal_drive` evaluates
an explicitly authored pressure, head or injection law inside the existing
native regional thermal residual. The pure callback receives the native
candidate thermal endpoint, the strict interval midpoint as its time coordinate,
and the original radial graph re-prepared at each pipe's actual candidate
temperature. The same original material law, support and radial history remain
in force; the adapter supplies no surrogate resistance or inferred pump law.

Every prepared residual, numerical derivative probe and final acceptance query
resolves its own drive. The original midpoint momentum step restarts the same
accepted radial history under that trial's forcing. Its actual signed flow
supplies the existing enthalpy transport, and only its resolved viscous heat
enters the caloric source. Pressure and recovering head work remain separate
mechanical exchange. The accepted drive sample travels through the original
thermal payload with the accepted flow and heat, rather than being reconstructed
or queried again after convergence.

The prescribed-drive entry still queries once per native interval and reuses
that sample. Its `thermal_time_s` is None; feedback samples explicitly record
Some(end_s), while `sample_time_s` remains the midpoint time. Both entries use
the same full-versus-two-half adaptive controller and authored event splitting.
Only accepted fine-step samples receive physical accounting credit.
`drive_evaluations` counts all actual queries, including discarded comparisons
and terminal source refusals; existing native-attempt, thermal-endpoint and
branch-solve budgets bound the work.

A feedback-source refusal remains a terminal `Step/Resolve/Drive` error.
Original material preparation failures remain momentum-source errors, and heat
exchange errors retain their separate variant. Source failure, nonlinear failure,
work exhaustion and cancellation preserve the complete input window, including
after a locally accepted prefix. Scratch workspaces remain reusable.

This is first-order endpoint-temperature/midpoint-momentum coupling with fixed
reference density, fixed regional masses and zero hydraulic compliance.
Temperature-controlled pressure is an external mechanical law, not a new fluid
equation of state or thermal expansion model. No flow-dependent pump-curve solve,
shaft/controller history, moving valve topology, cavitation or changed free-node
injection constraint is supplied.

Four focused regressions cover prescribed-drive parity and query semantics;
original-viscosity head feedback against independent annular momentum and
nonlinear caloric balances, including recovering head work and nonbinary events;
late source/cancellation refusals, endpoint budgets and retry; and exact authored
event arrival under decimal interval caps in both drive modes.

The unchanged bounded original-source recipe compiled the network implementation
and passed all 62 selected tests, including these four regressions. It links the
original rheology, radial momentum, caloric transport, nonlinear solver and
numerical dependencies with an existing context-only cancellation facade.
Full-workspace Cargo and native asupersync/runtime integration remain unverified.
Pinned rustfmt parsed and formatted the four network files, and its subsequent
check passed; only whitespace and optional trailing commas changed, leaving the
physical source, tests and independent numerical oracles unchanged.


### Adaptive conjugate finite-cell convection

Complete source-conjugate heat assemblies can advance through the original
adaptive implicit heat controller. Original wall conductivity, eliminated
surface temperature, source-pinned fluid properties and reciprocal conductance
are refreshed in every native coarse/fine/retried solve. The accepted fine
trajectory retains the original finite stores and heat/entropy account;
rejected proposals contribute numerical work only. Native forcing events land
at their actual times, and final source samples follow target-time events.
The original source, correlation, quasi-steady layer and fixed-topology limits
remain active; temporal refinement does not establish empirical accuracy.


### Adaptive source-conjugate finite-gas convection

`ConjugateGasSurface::adaptive_isochoric_body_candidate` connects the original
conducting surface to the existing adaptive finite-gas wall-feedback window.
Each coarse, fine and retried proposal refreshes the actual surface and uses
its own duration for the material/film source-work audit. Original caloric
identity, pure-phase admission and conductivity support at actual initial and
final native body temperatures remain mandatory, including insulation. Only
accepted fine proposals contribute heat and retained source history. Existing
fixed-speed and prescribed-flow films supply the native wall-window callback
through their original isochoric coefficient API. Reference source identity,
bulk-density velocity, full-area and disabled body-node radiation rules remain
unchanged; no additional storage or independently integrated heat is added.

## Adaptive droplet evaporation throughout a connected finite gas network

The sourced droplet chamber transport adaptive module advances the complete
original gas map and its named droplet cloud over one positive physical duration.
Each trial calls the original transport-then-chamber owner: every declared
pressure, diffusive, porous, membrane, dusty-gas, layered dusty-gas and effusion
path samples its actual incoming gas state, then the original droplet films
sample the transported atmosphere. Phase-only and reciprocal sensible-heat
entry points share this controller and retain original source applicability.

One coarse interval is compared with two original fine intervals. The shared
material norm includes every gas node's temperature, total mass and union of
species inventories, plus each droplet's caloric state, mass and geometry.
An empty droplet population still evolves and refines the network. Only accepted
fine intervals contribute to the common material physical-window account; the
returned endpoint retains only the last fine source/transport receipts.

Original phase step-size limits, species withdrawal and typed nozzle/porous/
membrane equilibrium-crossing refusals may reduce the span. Other source,
caloric, conservation and cancellation errors stay terminal with their original
payload and attempted clock. Comparison, reduction and parcel/species ceilings
bound work. Explicit immutable trial maps and reusable source-row buffers avoid
hidden physical history; a failed request publishes no partial state.

Four regression groups compare exact manually composed fine steps, refinement
against 256 original intervals, remote-only transport with no droplets, and
bounded retry/cancellation/source refusal. Independent final gas/liquid mass,
species and energy checks accompany the physical scenarios. Compilation,
runtime tests and rustfmt could not run after the execution workspace disconnected.

This remains first-order split evolution with local error estimates, fixed
volumes, positive droplet seeds and prescribed relative speed/heat coefficients.
It adds no dry-out event, droplet dynamics, external reservoir, or coupled
entropy/passivity certificate.


## Multiple conducting bodies on one adaptive finite-gas window

`conjugate::multiple::adaptive_isochoric_bodies_candidate` binds each physical
body ID to its original `ConjugateGasSurface` and uses the original simultaneous
multiple-body gas controller. Canonical bounded ID lookup keeps distinct
material/caloric sources attached to their own bodies throughout coarse, fine
and retried proposals. Each source eliminates its own actual fluid-facing
temperature against the common gas endpoint. No sequential gas update or heat
redistribution is introduced.

The original layer-work audit uses each actual interval duration and one body
share of the prescribed source-work allowance. Initial and final states of
every proposal retain exact caloric-object identity, pure-phase and actual
conductivity support admission, including insulation. Full positive exposed
area or explicit insulation, disabled body-node radiation and original Biot
admission remain mandatory. The controller alone apportions prescribed heats
and retains accepted fine heat, entropy and original source reports.

These remain quasi-steady normal conducting paths with independent lumped
bodies and a well-mixed gas. There is no added layer storage, lateral spreading,
spatial gas field or empirical accuracy claim. Three G1/G3/G4 regression groups
cover independent nonlinear layer/film powers, simultaneous native half-step
agreement, body-order invariance, original source identity/domain refusal and
cancellation. Pinned formatting/parsing passed; native integration compilation
and runtime execution remain unverified.


## Original source films on quasi-steady straight fins

`finite_gas::fin::SourceFinGasSurface` composes the original sealed film owner
and original all-solid conductivity source with the constant-section fin BVP.
Fixed-speed and prescribed-mass-flow films retain their original property
conventions, source identities and bulk-density velocity. Uniform h is sampled
at the actual base/gas trial; k is held at the actual base value. The original
source admits the computed tip-to-base temperature interval, including every
interior conductivity knot. Explicit bounds limit neglected k variation and
the largest source-based transverse Biot number, using the minimum admitted
k across that interval rather than just its base value.

`StraightFinArray` distinguishes bare base, covered fin-root footprints and
actual wetted side/tip area. The original film must cover the wetted area.
The native owner receives one combined bare-base-plus-fin conductance divided
by the original base area; raw film hA is retained as evidence, not reapplied
as heat. Original gas inventory and declared base-body caloric/storage remain
authoritative. The fin material may differ from the base-body material.

The isochoric convenience path reuses the original nonlinear wall-feedback
transaction and endpoint source checks, with body-node radiation disabled.
The coefficient callback is also available to original shared-body/adaptive
owners under the same explicit geometry/source/model assumptions. These are
quasi-steady attachments with no inferred fin mass or heat capacity. No fin
array flow, transient fin field, spatially varying h/k solve, radiation,
contact resistance or mutual shielding is supplied.

Three source/native G0/G1/G3/G4 regression groups cover independent heat laws,
source retention in both directions and equilibrium, interior conductivity
extrema/domain limits, source-wide Biot refusal, cancellation/retry, actual
native energy exchange, area ownership and conductivity-sensitive cooling.
All changed sources passed the pinned parser; these cross-crate integration
tests remain uncompiled and unexecuted. The four standalone analytic fin-kernel
tests executed successfully as recorded in the fs-convection contract.

### Resolved transient Carreau networks with original flow-dependent head laws

`graph::rheology::transient::head` adds a source-retaining `FlowHeadSource` and `FlowHeadLaw`; `TransientCarreauGraph::step_with_head_laws` solves each original radial midpoint response together with its authored signed flow-dependent pressure head. Continuous nonincreasing head laws are admitted so the exact condensed conductance `G / (1 - G dH/dQ)` remains positive for the original graph pressure solve. A separate pressure residual tolerance and bounded source-query budget control the scalar head solve; every original radial trial consumes the existing branch-work budget.

The original source owns its inclusive signed flow interval, pressure, derivative and any reference-density restriction. No extrapolation, clipped source value, inferred reverse quadrant or substituted viscosity is permitted. Midpoint source queries remain inside support; both endpoints of the accepted affine flow history must also be supported. Final-flow support is checked after the coupled graph pressure/storage solve, so an unaccepted pressure guess cannot prematurely reject an otherwise admissible operating point. Mathematical head brackets do not require evaluating inadmissible extreme radial profiles.

Accepted reports retain the exact shared source identity, actual midpoint head and derivative, effective conductance, and signed `H(Q_mid) Q_mid dt` work. Numerical head-residual work remains separately visible in the mechanical reconciliation; it is never deposited as viscous heat. Existing held heads, radial inertia, original nonlinear rheology, compliance work, source errors and cancellation retain their native ownership and atomic state publication.

Verification: the bounded existing original-source network recipe compiled and passed 67 tests, including five new independent annular/nonlinear-head, Schur-tangent, extreme-support, coupled-storage-support and refusal/retry regressions. The recipe uses the actual production numerical modules with a context-only cancellation facade; full Cargo integration and native asupersync runtime execution remain unverified.

### Original pump sources in regional transient Carreau heat transport

The regional transient owner accepts original per-branch `FlowHeadSource` objects through `advance_with_head_laws` and the scheduled/thermal-feedback `RegionalCarreauDrive`. Each thermal trial restarts the same accepted radial history and solves the actual nonlinear equipment operating point against its original endpoint-prepared viscosity. The accepted drive and native head report retain the same shared source. Actual signed head work remains mechanical; only original integrated radial viscous heat is deposited. Source/material refusal, cancellation or a later trajectory failure returns neither replacement physical history.

Original head queries are bounded across every thermal residual and every coarse/fine comparison, including discarded trials. Completed work counts remain separate from accepted physical credits; a terminal refused native call is bounded by its remaining allowance without inventing a partial count. The original thermal inventory's `max_substeps` also caps the entire accepted trajectory across authored drive events.

Adaptive resolution additionally compares each branch's viscous heat, held-head work and original nonlinear-equipment work under explicit absolute/relative transfer tolerances. Signed transfers and sums of absolute half-step transfers are checked separately, preventing opposing actuator work from disappearing behind unchanged velocity or temperature. These are local estimates of the native accounted transfers, not a guarantee for unsampled forcing or separately cancelling mechanisms hidden inside a caller's heat payload.

`pump::radial::RadialPumpHead` retains the exact original `Arc<PumpHeadCurve>` and `Arc<ResolvedLiquidState>`. It delegates actual pressure and optional drive-power admission to the original equipment evaluator, and supplies its exact density/gravity-scaled piecewise-linear tangent (right derivative at an interior knot, left at the last endpoint). Head/drive support is intersected without extrapolation; the original nonnegative pump quadrant and viscosity limits remain enforced. Native constant reference density must equal the original inlet density. This is neither an expansion model nor an inferred reverse quadrant, affinity scaling, efficiency or motor-heat model.

Verification: the unchanged bounded original-source network recipe compiled and passed 72 tests, including independent nonlinear equipment/annular/caloric balances, signed recovering work, source identity, aggregate query refusal/retry, opposing quadratic drives at constant state, and original native substep exhaustion. Two older convergence fixtures now explicitly author sufficient native inventory allowance; their physical tolerances and oracles are unchanged. Pinned formatting checks passed for the changed modern modules and new adapter/test children, and focused diff checks passed. The numerical recipe uses actual production numerical sources with a context-only cancellation facade. The two original-material/pump adapter integration tests are authored and statically reviewed but outside that recipe; they and full Cargo/native runtime integration remain execution-unverified.


## Adaptive prescribed-flow convection on original finite inventory

`transient::nonlinear::adaptive` advances the original
`advance_implicit_with_heat` transaction over a bounded physical window.
Full-versus-two-half backward-Euler comparisons check every active cell's
specific storage coordinate [J/kg] and original T(h) [K], including
nonadvecting finite walls. Relative tolerances scale interval changes,
preserving formation-datum independence; enthalpy remains visible on latent
plateaus. Only original accepted half-step states, heat/source payloads and
audits are retained, with no extrapolated state or replacement caloric law.

The original pump/inlet schedule is sampled independently at each actual
endpoint and its knots are exact comparison boundaries. Resolver callbacks
receive the particular native interval's original start and candidate endpoint
for duration-dependent films. Original source/numerical/physical refusals are
terminal, including discarded coarse trials. Separate attempt, source-query
and accepted-history limits bound work; the original thermal max_substeps also
caps accepted steps. Compensated final per-cell storage is checked against
accepted external heat, with a separate sum of absolute native energy defects.
Thermal topology/boundary events remain explicit left-limit handoffs; fixed
mass/geometry, spatial refinement and forcing-resolution assumptions remain.

Seven focused native regression definitions cover scheduled analytic residence
response, latent cancellation, original half-step source equivalence,
nonadvecting wall cooling, work/source/event refusal and cancellation. They
were parsed with the pinned formatter but have not been compiled or executed
in this environment. This is a local temporal estimate, not a global error bound.


## Distributed nonlinear source fins on the original finite-gas owner

`finite_gas::fin::nonlinear::NonlinearSourceFinGasSurface` resolves the axial
fin field with the original all-solid conductivity curve and original gas-film
selection. Every edge consumes the curve's exact integral mean, including all
crossed source knots. Every side/tip node queries the original film at its actual
temperature and retains that typed source payload. Fixed-speed and prescribed
mass-flow constructors keep the original bulk-density velocity convention.
The initial uniform base field never requires solid conductivity at an unused
gas-reservoir temperature. Unsupported states refuse without extrapolation.

Bare base area, covered root footprints and wetted fin area remain distinct.
The accepted base side sample supplies the bare-base coefficient without an
extra source query; a single combined base-to-gas conductance enters the native
heat owner. Node film h*A values remain evidence only. Reports retain the
actual temperature field, original conductivity source, signed base/side/tip
heat and the array's signed and sum-absolute spatial heat defects. The fin
material may differ from the base body's material, and no fin mass is inferred.

The absolute heat tolerance covers the whole repeated array and is divided by
fin count before the representative spatial solve. The original finite-body
wall-feedback convenience path additionally caps that tolerance by its actual
interval source-energy allowance divided by duration, and disables relative
widening of that allowance. Native body-node radiation must remain disabled.
Initial, trial and actual accepted endpoints use the same source gates and
return the original body/gas transaction or a typed refusal, with no partial
physical publication. Cancellation retains the original checkpoint error.

This model remains one-dimensional and quasi-steady, with a uniform gas
inventory and explicit local film interpretation. The selected correlation's
authored characteristic length and physical domain remain authoritative;
local evaluation does not establish a resolved non-isothermal boundary layer,
fin-array interference, fan power, transient fin storage, contact resistance,
radiation or arbitrary nonlinear uniqueness. Transverse Biot admission applies
to evaluated local states, not an unsampled global constitutive bound. Spatial
refinement remains separate from nonlinear tolerance.

Three integration regression groups independently reconstruct original-source
face and surface heat, check source retention and support/refusal behavior,
and couple the profile to actual reciprocal body/gas energy. They were parsed
and reviewed but remain uncompiled and unexecuted in this environment. The
exact numerical fin module's ten tests, including six nonlinear groups, were
compiled and executed successfully; that is not a full fs-airflow/DSR pass.

### Original equipment power retained by native head evaluations

The existing `PumpPowerState` and `PumpLossRouting` value types now have one neutral owner in `graph::rheology::transient::head::power`; their established `pump::power` paths re-export the same types and preserve every field. `FlowHeadResponse::power` optionally retains the original evaluator's power state and its original loss partition. `RadialPumpHead` supplies that value from the same already-budgeted head/drive evaluation. Accepted heads therefore carry actual shaft/electrical input and separated losses without another source query.

Native admission checks each nonnegative pump/drive power channel against actual `H(Q)Q`, and each routed heat channel against its original fraction. Roundoff allowances scale only the same compared channel: large through-power or an unrelated loss cannot conceal invented small heat. Head-only and signed recovering sources remain valid with no power data; no efficiency or heat is inferred. Authored nonzero shutoff input/losses remain possible when hydraulic work and pipe dissipation are zero.

Verification: the actual standalone production power leaf compiled with the pinned compiler, and its large-throughput/small-heat regression passed (1 test). Original native graph regressions for accepted zero/nonzero-flow power, exact query counts and failed-partition rollback, plus original pump adapter receipt assertions, were authored and formatter-parsed but not executed after loss of the prior numerical dependency projection. Full Cargo/native runtime integration remains unverified. This native slice retains values; explicit regional heat destinations are supplied by the subsequent thermal composition.


## Finite fin storage on the original thermal timeline

`deposition::fin::TransientFinTopology` lowers one constant-section fin with
explicit axial partitions, masses, caloric identities and spatial placement
into the original native thermal cells and links. Root/internal contacts use
the physical cross section and explicit half-distances; lateral and optional
tip contacts own distinct surfaces. The root retains an authored base distance
and contact resistance. No base bare area, fluid inventory, density, heat
capacity, steady efficiency or second thermal integrator is inferred. Identity,
geometry, boundary, work-budget and cancellation refusals publish no topology.

Five focused regression groups compiled and executed successfully against the
exact topology and original native thermal/conduction routines in a bounded
source projection. They cover finite fin startup storage, independent
three-store eigenmode refinement, mass sensitivity, energy closure, nonuniform
geometry and coverage, insulation/authored sources, and transactional refusals.
The analytic temperature errors were 0.305197664745, 0.151894055408 and
0.075771342836 K for steps 0.02, 0.01 and 0.005 s. Dependency-heavy unexercised
phase and implicit module registrations were omitted; no exercised thermal
updates were replaced. This is not a full workspace, source-film integration
or DSR pass, nor experimental fin validation. Spatial refinement and transverse
uniformity remain independent physical obligations.

### Original equipment input and routed losses in regional Carreau heat evolution

`TransientRegionalCarreauTransport::advance_with_drive` borrows the complete original `RegionalCarreauDrive`. Its optional `branch_heat_destinations` maps each branch's pump and drive losses to explicit native `CellId` recipients through `EquipmentHeatDestinations`. Positive original routing fractions require recipients even when the current loss is zero. Missing, foreign or malformed destinations refuse; destination data without original accepted input-power data also refuse. Head-only sources remain usable without invented input power or heat.

The shared endpoint owner reads each accepted head's exact original power/routing receipt. It performs no second equipment query. It deposits only the authored fractions of pump and drive loss, retains the original surroundings fraction, and keeps the original radial viscous heat separate. Distinct recipients may represent different coolant regions; shared recipients accumulate both components and any independent existing heat. Zero hydraulic flow can therefore heat finite inventory from authored shutoff pump/motor consumption. The same native caloric residual then feeds local temperatures back into original viscosity and the equipment operating point.

The exchange callback continues to own additional heat and streams. Its source values and native ambient/boundary heat are measured before equipment losses are added, and remain `other_heat_j`; automatically routed losses must not be added again by the callback. No heat is inferred from hydraulic work, unknown efficiency or a rejected numerical residual.

Per-branch reports retain hydraulic, shaft and electrical work, both losses, separately routed pump/drive heat and surroundings heat. The accepted endpoint passes both the hydraulic-plus-deposited-heat balance and an independent boundary balance replacing only known equipment hydraulic work and losses with original electrical input minus surroundings heat. Held/head-only sources remain explicit mechanical inputs. Adaptive comparisons additionally resolve every original input/loss/routed channel independently; discarded coarse or rejected trials consume work but receive no physical credits. The whole requested window remains atomic, with the original source-query, radial-work and thermal-substep caps.

This composition retains fixed original density/masses and original caloric identity; it supplies neither expansion nor shaft dynamics. Time staging remains original endpoint thermal feedback with native midpoint momentum/source quadrature.

Verification: the pinned formatter parsed and checked all eight integration files. Four focused native regression definitions cover zero-flow heating with distinct/shared recipients and independent callback heat, nonlinear equipment/annular/caloric feedback against an independent scalar/matrix oracle, time-varying electrical input resolved while cooling holds all states fixed, and missing/foreign/late destination refusal with retry. These integration tests were statically reviewed but not compiled or executed after the previous numerical-source projection was lost. The separately compiled original native power leaf's small-heat/large-throughput regression passed; no full Cargo/native runtime or current regional numerical test pass is claimed.
