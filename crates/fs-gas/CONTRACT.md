## Acceleration through reacting and molecular network compositions

The native `ReactingNetwork1d` and `DissipatingNetwork1d` owners now
accept `NetworkBodyForces` through their original complete-duration retry
loops. Chemistry, diffusion, heat conduction, viscosity and finite-port
transport retain their original state and constitutive laws. Distributed
acceleration acts inside the transport interval at the actual common clock;
if it shortens the acoustic step, all enclosing chemical/molecular/wall/valve
halves restart from accepted fields and finite-plenum inventories.

The explicit molecular entries and original RKL2 entries share this same
transaction. RKL2's RHS work allowance spans both halves and discarded
attempts; body-force histories never use its artificial stage index as time.
The existing reservoirs, finite walls and valves remain available. The
reacting composition also retains prescribed plenum volume rates and their
separate pressure work. No extra source model, duplicate gas owner or
source/chemistry energy correction is introduced.

`transport.body_force_work_j` is removed exactly once from the complete
referenced and physical energy residuals. Chemical formation-reference
changes, external species enthalpy, plenum pressure work and wall heat/work
keep distinct accounts. The final accepted chemical/molecular state and
every source report publish together. All-zero body forcing returns to the
corresponding original numerical path. Source Mach admission applies at the
mechanical source endpoints/extrema; the original chemical/molecular
temperature and constitutive domains still govern those separate stages.
These APIs do not combine the explicit global entropy gates or adaptive
molecular error controller with nonzero body forcing. RKL2 extends diffusion
stability, not acoustic limits or a global accuracy/positivity certificate.

Five focused regressions cover force-history/chemical replay after a shortened
interval, physical energy in accelerated reacting molecular networks with
baths and compression, shared RKL2 work-limit/late-cancellation refusal,
explicit/RKL2 molecular source work and exact inactive-path replay after
active forcing. Native compilation/execution was unavailable during this
implementation; these regressions are not claimed passed.

## Distributed acceleration in the native multicomponent network

`MixtureNetwork1d::step_with_body_forces` and
`step_with_body_forces_and_boundaries` apply one original `BodyForce` law
per retained branch, including cellwise acceleration/force density and absolute-time
linear histories. Exact mechanical half maps surround the original simultaneous
SSP RK2 branch/plenum transport. The first accelerated field supplies every
characteristic port and the SSP convex average; the second source uses actual
transported mass. Source work is represented impulse paired with midpoint velocity,
with original species calorics, Mach admission and interior-history speed extrema.
Full-interval CFL and all retries belong to the common network clock. A shortened
trial restarts all branches and plenums from accepted storage and reevaluates the
history on the shorter physical intervals.

Accepted branch reports retain axial source impulse/work and representation
defects. The network sums energy only: `body_force_work_j` is subtracted once
in its first law, separately from exterior enthalpy and compact-volume pressure
work. `body_force_energy_residual_j` remains visible and is never heat or a
residual correction. Internal finite-port species/energy exchanges remain exact
opposites; unrelated branch axes are not summed into a fictitious global momentum.
Prescribed reservoirs and moving compact volumes use the existing boundary laws.
Plenums remain stationary perfectly mixed owners, and their supports retain port
impulses; no moving-frame plenum momentum model is implied.

Active forcing reuses each branch's original lazy source buffer and bounded
construction-time receipts. Zero forcing uses original transport arithmetic even
after those buffers have been used. State, volume, clocks and accepted reports
publish only at the existing common commit point. This is the estimated native
source/transport composition, not a well-balanced hydrostatic scheme or a
continuous peak certificate for general histories. The explicit entropy-admitted
network entry points remain separate until complete source endpoint entropy is
accounted. Focused regressions cover analytic periodic histories/work, actual
stage mass with finite ports, simultaneous reservoirs/compression, hidden CFL
peaks, zero-force replay and late cancellation. Native compilation/execution was
unavailable during this implementation; these regressions are not claimed passed.

## Physical-time acceleration and force-density histories

The original BodyForce owner now admits uniform and per-cell linear histories
for acceleration or force density. Each history specifies an absolute reference
time, its value there, and its constant time derivative. The shared kick
integrates the two actual source half-intervals separately: [t,t+h/2] and
[t+h/2,t+h]. It carries the accepted clock and within-step offset separately,
so a shortened retry or a later reacting/molecular window cannot reuse a stale
requested-interval sample. A source-only apply starts at the current flow clock
and does not advance that clock.

Frozen caloric, evolving-species, and retained NASA/Peng--Robinson gas owners
use the same impulse kernel. Acceleration uses each stage's actual mass density;
force per volume retains its independent volume weighting. Actual represented
momentum change supplies conjugate kinetic work, and representation loss stays
in the original residual. No chemical energy or heat is synthesized. Old and
moved cell volumes remain owned by the native ALE stages, including reciprocal
inertial-wall and nonideal piston solves. Time-varying force density on changing
volume retains the native second-order source/transport composition; exact
integration of the source history alone is not an exact moving-volume or
coupled-flow solution.

An affine force can reverse inside a source stage even when its net impulse
vanishes. The shared kernel locates that strict interior velocity extremum;
each original EOS decodes the represented state and applies the authored Mach
bound. Stationary and moving relative-wave CFL gates also admit the peak.
A piston source before the wall root bounds every speed in the original
declared bracket. Interior probes add no impulse, work, state buffer, or source
receipt. Ordinary held sources and zero-rate histories retain their previous
arithmetic; zero histories still use the unforced path. Refusal and cancellation
preserve accepted gas, geometry, wall state, source accounting, and clock.

# fs-gas contract

## Conservative finite-plenum exchange candidates

The existing frozen and evolving-composition network `Plenum` types now expose
`exchange_candidate` for composing their finite storage with other native gas
owners. Positive mass/species masses [kg] and total energy [J] enter a rigid
volume. The update advances extensive inventories, then recovers density and
internal-energy density at that same volume. Both endpoints pass the exact
original `FrozenGas` / `MixtureGas` caloric admission. Pressure and temperature
therefore respond to filling, blowdown and transported energy; no pressure
source, fixed gamma or second EOS is introduced.

The energy transfer is the port's complete referenced total-energy exchange.
It already includes flow pressure work. Incoming kinetic energy becomes
internal storage because compact bulk momentum remains zero; the composing
owner must account for the opposite Cartesian/axial impulse at its support.
Species retain canonical order and their original formation-energy gauges.
There is no extra enthalpy, reaction, heat or moving-volume-work term.

Frozen plenums are `Copy` stage values and allocate no storage.
Mixture `exchange_candidate_into` prepares a bounded stack candidate and
publishes into an existing species vector only after the final cancellation
poll; successful repeated calls allocate no stepping storage. The convenience
mixture return-value method allocates one bounded species vector. Exact zero
transfers preserve the original state. Species may reach exact zero, but
exhausted total mass, negative inventories, unsupported caloric states,
nonfinite exchanges and unrepresentable extensive storage refuse without
clipping or partial output.

These are shared algebraic stage operations, not independent time integrators.
A coupled owner supplies its common-stage fluxes, compact-volume CFL, retries
and joint publication. The existing 1-D network integration arithmetic is
unchanged. This does not resolve internal tank waves, mixing gradients, neck
inertia or molecular/thermal interfaces, and does not duplicate the distinct
material-parcel/NSF-DG coupling.

Five focused native tests passed using these actual production leaves,
the original Plenum definitions and original EOS/characteristic dependency
closure. They cover variable-Cp filling and incoming kinetic energy, large
formation-gauge invariance, source-port blowdown cooling and pressure feedback
with independent adiabatic/RK2 refinement checks, canonical species mixing and
depletion, and transactional source/shape/exhaustion/late-cancellation refusal.
Unrelated network/thermal/chemistry wrappers were outside that bounded source
build; this is not a full-crate or workspace Cargo/DSR gate.

## Compact stagnation reservoirs on planar and volumetric gas traces

Frozen and evolving-mixture gas owners expose conserved and reconstructed
planar/three-dimensional compact-reservoir fluxes. They reuse the original
one-sided reservoir characteristic and source-owned choking law. Reservoir
pressure and temperature are stagnation data for inflow and pressure is the
outflow backpressure. A held exterior HLLC state remains a different model.

Incoming gas takes the stationary reservoir's original composition and
referenced stagnation enthalpy, with zero transverse velocity. Departing gas
retains interior species, both shear velocities and complete kinetic energy.
Pressure work is included once in transported enthalpy. Source-domain and full
caloric-precision admission remain binding; mixture output publishes only
after final cancellation. Successful face evaluation uses bounded stack
storage and performs no heap allocation or replacement caloric construction.

Five new native regression groups check exact gamma-law sonic inflow,
variable-Cp reversal and discharge, both orientations, original energy gauges,
different donor compositions and absent species, shear/energy transport,
conserved/reconstructed parity, and late cancellation with unchanged output.
All 18 selected source-linked native boundary tests passed. This bounded build
uses the actual boundary sources and original dependency artifacts; unrelated
network, thermal, chemistry and entropy-flux wrappers are excluded. It is not
a full crate/workspace Cargo or DSR result.

This is a stationary compact reservoir with externally held state. It supplies
neither finite mass/energy feedback, resolved mixing, entrance loss, neck
inertia nor molecular reservoir data. A future finite-volume owner must evolve
the reservoir with the exact opposite area-integrated port transfers.

## Planar and volumetric acoustic pressure outlets

The frozen and evolving-mixture gas owners expose characteristic outlet fluxes
for both conserved and reconstructed planar/three-dimensional traces. They use
the same original-caloric incoming-characteristic solve as the 1-D outlet,
with reference pressure and reference velocity signed along the supplied normal
axis. The outgoing entropy, composition and each shear velocity are retained;
species, momentum, pressure work and every kinetic-energy component cross the
same face. Shared outgoing-flux assembly avoids a second constitutive closure
and allocates no heap workspace. Full caloric precision, source-domain limits,
sonic sampling, supersonic causality, backflow refusal and cancellation remain
binding. Mixture output publishes atomically, including late cancellation.

Focused G1/G3 regressions check exact outgoing simple waves, dimensional
reduction, both normal orientations, variable source Cp, absent species,
both shear energies and sonic/supersonic regimes. These are local smooth
normal-wave conditions; they do not guarantee nonreflection for oblique waves,
shocks or interactions with transverse composition/entropy gradients.

## Adiabatic throttling with original phase enthalpy

`material::phase_throttle::throttle_adiabatic_mixture_stream` closes the steady
pressure-reduction process with zero original enthalpy increment and no shaft
work or external heat. Its feed is the same original-source/temperature/
pressure/mass-flow type accepted by the mixture heat-exchanger and mixing
owners. The existing PH source resolves Joule--Thomson temperature change and
phase redistribution with retained species data and phase policy. Formation
and entropy reference offsets cancel before rounding; no constant-Cp,
isentropic-expansion or separate latent-heat approximation is introduced.

The supplied mass flow is conserved. Actual outgoing-minus-incoming enthalpy
flow and entropy production are checked in W and W/K before publication;
stream pressure work is already included in enthalpy. The source energy
tolerance is bounded by the physical power allowance divided by throughput.
The existing source work limit includes the inlet selection and every outlet
trial, with no second accepted-state evaluation. Equal pressure has exact
zero energy/entropy change after one original PT query. Cancellation, invalid
pressure direction and source/bracket/work refusals publish no partial outlet.

Focused `phase_equilibrium_throttle` regressions cover the independent
cryogenic PR phase/enthalpy reference, the ideal-gas dilute limit, pressure-flow
work, reference-datum and throughput invariance, exact rest and refusal/replay.
These verify the declared model, not a measured valve or refrigerant dataset.
The process assumes negligible inlet/outlet kinetic and potential energy
changes, steady flow and phase equilibrium. Mass flow and downstream pressure
are inputs; valve geometry, discharge/choking, slip and phase kinetics are not
predicted. The thermodynamic convention follows
[MIT 5.60, lecture 4](https://ocw.mit.edu/courses/5-60-thermodynamics-kinetics-spring-2008/e1f8054833e8077bf9428c27cbf0465c_5_60_lecture4.pdf).

## Original real-gas pressure at impermeable end walls

With `material,experimental-euler`,
`RealGasBoundary::PhysicalReflecting(RealGasWallControl)` and
`RealGas::reflecting_wall_flux` resolve stationary wall traction using the
retained NASA/Peng–Robinson source. Positive velocity toward a wall uses the
existing original-law piston Rankine–Hugoniot solve and its conservation,
entropy and sonic admission. Retreating gas uses the original finite-entropy
inversion along `s = constant` and solves
`|u| = integral(rho_wall..rho_inside) a(rho,s) d(rho)/rho`.
The logarithmic density coordinate removes the singular-looking `1/rho`
factor. This produces a positive source pressure for strong expansions where
the legacy mirrored HLLC/HLL traction can become tensile.

Both sides return exactly zero mass and energy flux. Pressure enters the same
face momentum balance and wave-speed CFL accounting as the interior fluxes.
The physical boundary is evaluated at both SSP-RK2 stages; native field retries
and cancellation still publish no partial accepted state. Left/right wall
orientation changes the signed incoming velocity, not the positive traction.
At equilibrium the pressure is exactly the original EOS pressure without an
additional wall solve. Nonzero waves within the authored absolute tolerance
or the floating-point weak-shock resolution use its acoustic impedance
tangent. Every such continuation still admits the native isentropic density
and temperature path, including strict endpoint direction/bracket checks. The
authored pressure tolerance gates acoustic consistency with the native
isentrope at the linearized density; it is not a certified error bound against
the fully resolved wall wave. Small velocity does not bypass a source
temperature or density boundary. The compression search extends its
lower ratio toward unity as needed, and tightens caloric/pressure residuals
to resolve the quadratic sonic separation; a fixed lower bracket therefore
does not leave an unsupported gap next to the acoustic limit.

Controls retain explicit compression and rarefaction source paths, caloric
temperature bounds, velocity and quadrature tolerances, and finite work limits.
The rarefaction endpoint must itself be on a supported homogeneous isentrope,
even if a shorter wave would fit. The embedded three/five-point Gaussian
quadrature uses a fixed stack, bounded depth and a divided error budget;
original material evaluations retain their own allocations. Energy/entropy
datums cancel in the original finite-increment owners. Source, bracket,
quadrature or root exhaustion returns a refusal rather than pressure clipping,
an ideal-gas fit or an extrapolated vacuum state.

This is an Estimated classical single-phase shock/rarefaction boundary, not an
entropy-stability certificate or a treatment of phase transitions/nonclassical
waves. The caller must declare a classical rarefaction path with monotone
characteristic speeds for the endpoint CFL bound to cover the entire fan;
native material stability does not certify that assumption or detect every
nonclassical source. Existing `Reflecting` explicitly retains its numerical
mirrored-flux semantics. Physical wall controls must be selected to cover the evolving
boundary-cell states; the native source still checks every evaluated state.
Focused regressions in `material_nonideal_wall` cover an independent dense PR
shock, an independently integrated dense rarefaction with a tensile numerical
mirror, orientation, original datum invariance, the acoustic limit, closed
field pressure/mass/energy balance, acoustic source-endpoint refusal, a sweep
across both weak-wave transitions and bounded-path/cancellation refusal.


## Original real-gas heat exchange in compressive fields

With `material,experimental-euler`, `RealGasEuler1d::step_with_walls` advances
the original homogeneous Peng–Robinson field and then applies the shared
implicit wall kernel over the actual accepted wave duration. A prepared
`RealGasWallOperator` binds one contact per original cell, the unchanged
material/context/energy gauge and physical mesh. Each local source freezes
its post-transport density and uses the original finite internal-energy
increment and actual density-dependent Cv, including PR departure calorics.
Neither constant Cv nor ideal-gas Cp-R replaces that original derivative.

Reservoirs, authored constant-capacity walls and retained original finite
Cp/h/u laws share the existing bounded heat solve; latent inventories remain
energy coordinates. Source-bearing finite walls retain their declared storage
and work convention. The gas total-energy account, wall h/u storage, external
reservoir heat and prescribed sliding-wall work remain separate. Relative-slip
drag heating is retained once in gas internal energy and adds no external heat.
The whole-field storage residual subtracts actual boundary transport and drive
work; it does not reinterpret wall enthalpy as rigid-body internal energy.

Wave and source candidates share the original retry/commit boundary. Gas
states, pressure/temperature/acoustic caches, finite-wall inventories, reports
and clock publish together only after complete source admission and the final
cancellation poll. Source or gauge/mesh mismatch, missing original wall law,
unsupported thermal endpoints and exhausted heat work refuse atomically. Wall
candidate arrays allocate at construction; no per-step field clone or second
thermodynamic solver is introduced. Splitting is first order, conductance and
drag remain authored engineering inputs, and this does not add spatial thermal
diffusion, evolving composition, phase flash or moving end-wall geometry.

The `material_nonideal_flow::thermal` native regressions specify independent
closed-form PR energy checks for heating/cooling and post-compression density,
nonlinear finite Cp and latent storage, datum invariance, sliding work and
late-cancellation/domain rollback. This extension was checked with the pinned
Rust parser and scoped source/diff review; these native regressions were not
executed because the required sibling workspace dependencies and remote build
tools were unavailable in this session.

## Energy-selective kinetic alpha escape

`AlphaEscapeTable` retains an identified, ordered table of nonnegative escape
hazards against physical energy per alpha. `KineticAlphaTransport::with_escape`
opts into the local loss law `dn_i/dt = -lambda(E_i) n_i`; absent or zero hazards
retain the closed population. Linear interpolation stays inside the original
table support. Every represented physical energy must be covered, including
accepted endpoints and volume-history knots. A moving energy coordinate uses
`E_i = E_i,initial (V_initial/V)^(2/3)` when sampling the law, not the original
reference-grid energy. No extrapolation, clipping or fitted rate is inferred.

The existing kinetic midpoint and step-doubling owners compose exact positive
exponential escape symmetrically with their collision and burn maps. Survivors
and small losses are evaluated separately; fixed-size stage arrays bound work
without allocating histories or rebuilding published grids. Compensated local
alpha-number, alpha-energy and companion-electron-energy increments enter
time-error acceptance independently. Only accepted fine trajectories publish
their inventories and accounts; domain/work refusal and cancellation are atomic.

The coupled frozen interval also limits the maximum represented escape exposure
to one e-fold whenever kinetic alphas exist or an alpha-producing channel is
enabled. Its maximum hazard is prepared once per source evaluation. Unresolved
steps return the original retryable work refusal, so both saturated coarse/fine
maps cannot falsely accept an empty stock under continuous birth or a spurious
collision/escape partition. The original attempt/minimum-step limits remain
binding: very stiff coupled laws may require smaller steps or refuse. The exact
standalone escape kernel and an empty source-off population need no such cap.
This guard does not turn the local error estimate into a global certificate or
claim stiff-uniform accuracy for the split system.

Escaped alphas carry their actual physical-bin kinetic energy out of the parcel.
They are neither thermalized in the plasma nor removed from Maxwellian ash.
Each also removes two companion electrons at the current electron pool's mean
sensible energy. This is an explicit quasineutral inventory closure: it retains
the admitted charge discrepancy, but supplies no ambipolar/sheath potential,
return-current work, momentum flux, pitch angle, spatial orbit or receiving wall.
Thermal helium, fuel supply/exhaust and thermal energy confinement retain their
separate original owners; confinement energy must exclude these carried exports.

`KineticAlphaBalance` exposes escaped alpha amount and kinetic energy, companion
electron amount, and its thermal energy. Its kinetic account becomes
`delta(Ualpha)+Halpha_h+Halpha_e+Ealpha_out-Ealpha_born-Walpha_on=0`.
The complete network account includes both energy exports and adjusts He4,
electron, baryon and charge inventories. These quantities are available to an
explicit receiving-body calculation; a plasma energy debit alone does not claim
wall deposition or a spatial loss spectrum.

Energy-dependent alpha confinement in a Fokker-Planck description is discussed
by Katsurai and Maekawa, *Nuclear Fusion* 17 (1977),
<https://doi.org/10.1088/0029-5515/17/5/005>. The table here is caller-authored;
it does not implement that paper's fitted assumptions or a magnetic transport
model. Tests cover exact selective losses/deep survivors, constant-hazard number
decay during collisional slowing, an independent compressing `lambda(E)=aE`
analytic limit, all-channel fueling/thermal-exhaust competition, restart and
atomic refusal, including a saturated continuous-birth case whose independent
positive equilibrium exceeds the requested amount tolerance. Numerical time
error remains distinct from energy-grid error
and the physical validity of the authored escape law. The reflecting collision
boundaries and absence of a kinetic-to-thermal tail handoff remain unchanged.

## Source-owned homogeneous real-gas fields

With both `material` and `experimental-euler`, `dynamics::nonideal::RealGas` binds the original `ThermochemicalPengRobinsonModel`, frozen mixture/context and one immutable reference density/temperature to the native `Conserved` state. Native density-energy inversion supplies pressure, temperature and acoustic speed. Internal-energy changes use the common original reference, including density-dependent residual energy; arbitrary source formation-energy and entropy datums do not change field dynamics. Referenced internal energy may be signed and is not an ideal-gas positivity test.

The existing const-generic HLLC/HLL kernel now accepts a constitutive admission callback. Both HLLC stars and each HLL fallback middle state are resolved through the original real-gas EOS. A source-invalid star can use the original bounded HLL speed expansion; cancellation propagates. The ideal/frozen-gas paths retain the same shared numerical algebra and admission semantics.

`RealGasEuler1d` reuses the original quasi-one-dimensional mesh, conservative volume update, compensated totals, step controls and balance reports. It admits smoothly varying face areas and producer-integrated physical volumes, nonuniform cell widths, periodic interfaces with one shared physical end area, zero-gradient traces, fixed reflecting walls and prescribed exterior states in the same energy gauge. Reflecting walls use the shared EOS-checked flux and exactly zero mass/energy transfer. Numerical scratch is preallocated and accepted primitives are cached. Every first Euler stage, second Euler endpoint and final SSP-RK2 average is source-admitted before the final transactional swap; actual face signal speeds determine both stage CFL gates and bounded time-step retries.

Spatial reconstruction accepts explicit `FirstOrder` or source-admitted `Muscl` through the existing `StepControl`. MUSCL limits density, velocity and temperature on physical nonuniform cell widths, then constructs BOTH traces with the original NASA/Peng--Robinson state owner and common energy gauge. One scale contracts the complete primitive increment when either original-source endpoint refuses; twelve contractions bound the work before the exact accepted cell is reused. No individual thermodynamic coordinate is clipped. Constant and unrepresented trace increments reuse the original conservative payload and primitive cache without constitutive reevaluation.

Two traces per cell are cached in one lazy reusable buffer, so each face consumes already-admitted states without a duplicate material solve. Nonperiodic end cells retain their original first-order boundary relation, including physical reflecting walls. Periodic fluxes remain single evaluated and shared; nonuniform geometry uses adjacent physical widths rather than subtracting distant center coordinates. The native conservative face update and SSP-RK2 averages remain unchanged. Accepted reports count source-limited cells and declare any whole-step first-order retry. A CFL-only retry keeps MUSCL; other source/stage refusals halve the interval and allow the more dissipative first-order path. Gas, source work, wall accounts, primitive cache and time still publish together.

The method is second order on smooth admitted fields; extrema, source limiting, end cells and first-order retries limit that statement. The caller declares one homogeneous stable density/temperature branch, and every evaluated state is checked. No phase flash, global entropy or positivity theorem, spatial viscous/thermal transport or evolving composition is inferred. Prescribed constant-area end motion has the separate moving-wall entry described below. The original `material_nonideal_flow` regressions retain caloric/acoustic source retention, datum-invariant compression, nonuniform free streams, closed-box balances, an independent dense shock and atomic refusal. The added `material_nonideal_flow/reconstruction` cases cover dense acoustic momentum refinement, first-order diffusion comparison, original-reference invariance, uniform nonuniform-mesh flow, pressure impulse and cancellation. Eight native `reconstruction/linear` tests exercise nonuniform affine exactness, monotonicity, coupled source admission, work bounds, cancellation, exact original-payload reuse and smooth conservative transport refinement.

## Fully-ionized fueling and selective thermal ash exhaust

`FusionParticleExchange` supplies explicit constant ion rates in mol/s and
first-order thermal-particle removal hazards in 1/s. Every array uses the
physical order `p, D, T, He3, He4`, independent of canonical material ordering.
`DtFusionTransport::particles` selects this optional source in the original
D-T, four-channel network and kinetic-alpha entrypoints, including prescribed
volume histories. The D-T bank refuses nonzero proton or He3 ports. No source
is inferred when the option is absent; the original closed-particle path stays
available.

Each ion brings or removes its own charge's companion electrons. Injected
ions and electrons have separate authored temperatures, with sensible internal
energy evaluated from their original monatomic source cards once per solve.
Incoming source coverage and the declared nonrelativistic temperature limit
remain mandatory, even when a tiny injected population would dilute into an
admitted bulk state. Outgoing particles carry the actual mean thermal energy
of their pool. Thermal fuel dilution, charge balance, changing heat capacities,
reaction competition, radiation and alpha slowing therefore see the evolving
populations. The supplied ionization state is assumed; neutral/pellet ablation,
ionization energy, bulk momentum, ambipolar potential, recycling and reservoir
dynamics are not supplied. These are local inventory sources, not flow-face
enthalpy fluxes; no additional inlet/outlet pressure work is inferred.

The common bounded source kernel composes exact exponential exhaust with
constant injection symmetrically. Positive survivors and small losses are
evaluated separately, retaining deep trace inventories and tiny transfers.
The original midpoint/step-doubling owners resolve simultaneous supply,
selective removal, burn and thermal feedback. Trials allocate no source history
and retain independent compensated particle/energy accounts; local transfer
error includes the compensated low parts so a weak late flow is not hidden by
a large accumulated account. Every candidate retains original material and
physical-domain checks; refusal or cancellation publishes no partial state.

Kinetic He4 is explicitly excluded from this thermal exhaust. The Maxwellian
He4 inventory now evolves as its own retained stock through supply, exhaust
and restart. Fast alpha number, distribution energy and pressure remain under
their existing birth/collision owners. The reflecting lower-energy alpha grid
has no thermalization handoff: even its low-energy tail remains kinetic and is
not silently reclassified as removable Maxwellian ash. Incoming He4 enters
the thermal stock and does not seed the kinetic distribution.

Both result balance families expose `particles: FusionParticleBalance`, with
all five ion transfers, companion electrons and separate incoming/outgoing
heavy/electron energies. Species, baryon and charge audits include these
external transfers. The full energy account adds particle energy out and
subtracts particle energy in. `FusionEnergyConfinement` remains a separate
non-particle energy term and must exclude energy already carried by these
ports. The published transfers allow an explicitly coupled receiver to debit
or credit its own finite inventory without inventing that receiver here.

A particle substep resolves at most one removal e-fold for each populated
thermal species or nonzero incoming source when supply, nuclear reactions or
another thermal source also acts. Otherwise saturated coarse/fine exhaust maps can
both erase the sustained solution `N=S/lambda` while reporting a small local
difference. The physical drive is retained across both particle substeps;
reaction products created between them are checked before the second exhaust.
Pure unforced exhaust keeps its unrestricted exact map. Empty sinks do not
impose this restriction, and kinetic He4 is never counted as removable thermal
ash. The original attempt and minimum-step budgets remain binding; insufficient
resolution returns `StepBudget` without publishing an accepted result.

Fuel replenishment and helium removal are distinct physical needs in a
sustained burn; see the [UKAEA PROCESS fuel and ash exhaust discussion](https://ukaea.github.io/PROCESS/eng-models/cryostat-and-vacuum-system/).
The rates here are caller-authored closures, not a fitted pumping or magnetic
particle-confinement model. Numerical checks cover analytic supply mixing and
deep exhaust, initially absent fuel and restart, all five ion ports against an
independent continuum energy ODE under compression, second-order refinement,
all-channel Coulomb/radiation/confinement coupling, thermal-only kinetic ash
removal, cold helium supply without artificial fast alphas and atomic refusal.

## Finite thermal energy confinement and neutron-source handoff

`FusionEnergyConfinement` explicitly authors separate ion/electron energy
relaxation times and exterior temperatures. `DtFusionTransport::energy_confinement`
applies that closure through the common thermal kernels used by the original
D-T burn, the four-channel network, kinetic-alpha evolution and prescribed
volume histories. For fixed populations, each pool obeys
`dU_s/dt = -(U_s - C_s T_exterior,s)/tau_s`. Exact exponential source maps
retain the positive survivor under stiff cooling and use `expm1` for small
transfers. Symmetric composition retains the existing local error controls.

Separate signed `heavy_transport_energy_out_j` and
`electron_transport_energy_out_j` enter both error estimation and independent
publication balances. Positive values leave the plasma; a hotter exterior
produces negative outgoing heat. These are distinct from neutron/photon escape,
Coulomb redistribution, imposed powers, alpha storage and mechanical work.
The full account adds both transport exports to the existing energy balance.
Kinetic alphas lose energy only through their original collision owner; this
closure removes no particles and supplies no fast-alpha escape model.

The conventional zero-dimensional `P_transport = W/tau_E` interpretation and
separate radiation accounting are described by
[UKAEA PROCESS](https://ukaea.github.io/PROCESS/physics-models/plasma_confinement/).
Separate pool times and finite exterior-temperature offsets are explicitly
authored extensions here, not a fitted magnetic confinement or turbulence law.
The exterior is prescribed; its own thermal storage is not evolved by this
closure. Numerical tests cover analytic relaxation, signed reservoir heating,
stiff and tiny exposures, burn feedback, restart, compression/work refinement,
kinetic-alpha accounting and atomic refusal. `compressed_fusion` compares the
existing rigid and compression cases with an explicit finite-confinement case.

When a pool also receives imposed power, nuclear/alpha heat, reciprocal
exchange, radiation or particle exchange, each frozen source interval resolves
at most one of that pool's confinement e-folds. Saturated full/two-half cooling
maps can otherwise agree after erasing a continuously supplied thermal store.
For constant power and an initially equilibrated pool, the exact solution is
`U-Uexterior = P*tau*(1-exp(-t/tau))`. Unforced relaxation retains its
unrestricted exact map, and an inactive pool's shorter time does not constrain
an independently heated pool. The original attempt and minimum-step budgets
apply; unresolved coupled sources return `StepBudget` without publishing
state. This guard prevents the false saturated estimate; it does not replace
temporal refinement or establish a global error bound.

All three result families expose `neutron_emission()`, returning the shared
L1 `FusionNeutronEmissionV1`. This preserves the original accepted D-T and
D-D-neutron extents, product kinetic energies and interval duration. It leaves
the plasma's neutron energy debit unchanged. Downstream owners choose their
emission-time and energy projection explicitly; no deposit or spectrum is
inferred from an escaped-energy scalar.

## Accepted-interval fusion neutron histories

`DtFusionSolve::record_neutron_history` (also shared by `FusionNetworkSolve`)
opts into one `FusionNeutronIntervalV1` per accepted adaptive interval. It is
false by default. The D-T, network and kinetic-alpha results expose the records
through `neutron_history()`: disabled recording is `None`, while a recorded
zero-duration solve is `Some` with an empty slice. Positive-duration intervals
with zero emission remain present, so the accepted clock has no artificial gaps.
Both rigid and prescribed-volume entrypoints use the same option; moving-volume
records end at every original geometry knot.

Each trial keeps separate compensated D-T and D-D-neutron amounts and kinetic
energies, credited directly by the original finite-stock reaction maps. The
history never differences cumulative source accounts: a representable weak late
source survives even when it is smaller than an ULP of the previous total.
Midpoint predictors, coarse trials and rejected attempts are discarded. Only
the accepted two-half solution contributes a record. Kinetic fusion obtains
these moments from the same shared reaction owner while alpha births, collision
heat, storage and pressure work retain their existing independent accounts.

Record boundaries are the native represented relative-time endpoints; their
emission duration is exactly `end-start`. Separate interval reductions reproduce
the cumulative accepted source to roundoff, without requiring bit equality of
differently grouped sums. Recording does not alter source evaluation, local
error control, accepted step selection or the physical state. Records supply
integrated number/energy moments, not an intra-interval waveform or neutron
spectrum; a downstream time/energy projection remains explicit, with its own
resolution and domain requirements.

Numerical trials remain allocation-free. Only requested accepted records may
grow a vector, with fallible allocation and at most `maximum_steps` entries.
The optional history stays private until final source/state publication and
cancellation checks succeed; any refusal returns neither a partial state nor a
partial history. Focused regressions cover changing D-T source intensity,
rejected-step exclusion, late cancellation, empty histories and zero-emission
intervals, both neutron channels, geometry boundaries, a late source below
cumulative rounding resolution, and unchanged kinetic/confinement evolution.

## Kinetic fusion through prescribed compression and expansion

`fusion::kinetic::advance_adiabatic_kinetic_fusion_volume` moves an original
hybrid state without advancing its clock or redistributing alpha populations.
For the explicitly confined, co-moving isotropic nonrelativistic closure,
momentum scales as `V^(-1/3)` and every physical alpha energy as `V^(-2/3)`.
The shared stable adiabatic map returns separate work on thermal ions, electrons
and kinetic alphas, with positive sign for compression. Source temperature-path
admission includes absent products and gaps. Original source identities, fixed
thermal He4 stock and independent kinetic inventory survive the transformation;
the kinetic transport supplies its explicit nonrelativistic and dilute limits.

`advance_kinetic_fusion_with_volume` couples the existing `FusionVolumeHistory`
to all four light-ion channels, finite-time alpha drift/diffusion, background
Coulomb exchange, radiation and held external powers. During trials the alpha
grid has one scale `s=(V_initial/V_stage)^(2/3)`. Physical centers and midpoint
cell widths are evaluated in bounded stack storage; geometry changes no bin
amounts and performs no remap or artificial distribution broadening. Each stage
uses its actual volume for density, its physical grid for collision coefficients
and heat rewards, and its source-event birth energy for conservative injection.
Only publication allocates the final physical grid. Restart reanchors that
physical grid at its new volume rather than discarding prior compression.

Adiabatic maps surround the midpoint source solve, and attempted intervals stop
at every authored knot. Thermal and alpha work increments enter time-error
control. Accepted endpoints, including expansion/compression turnarounds,
admit complete source coefficients and birth support without rebuilding a
collision matrix. Thus supported neighboring midpoints cannot hide a knot whose
expanded grid no longer covers an enabled birth. Every physical stage also
checks the kinetic energy ceiling and the separate radiation motion limit.
Unsupported grids or sources refuse without clipping or extrapolation.

The alpha balance is `delta(Ualpha)+Halpha_h+Halpha_e-Ealpha_born-Walpha_on=0`.
The full account is
`delta(Uh+Ue+Ualpha)+En+Ephotons-Q-Eexternal-Wh_on-We_on-Walpha_on=0`.
The result retains the original history, three work accounts and total work.
Returning to the initial volume need not give zero net work when reactions or
heat transfer change the populations and energy during the cycle. The rigid
entrypoint reports zero work and agrees exactly with a constant two-knot history.

The momentum/energy law follows Parker (1965), printed p. 15, Eq. (2),
<https://polaris.nipr.ac.jp/~ryuho/coem/Parker1965PSS.pdf>; the distribution
need not be Maxwellian when isotropy is retained (Wierzchucka et al., section
2.1, Eqs. 2.19–2.22, <https://arxiv.org/abs/2603.25669>). Boundary motion,
confinement and isotropy remain prescribed assumptions. This scalar geometry
does not supply piston inertia, shocks, anisotropic magnetic compression or
spatial heat/particle transport. Time-error control remains separate from
finite energy-grid and model error. The `compressed_kinetic_fusion` example
compares rigid and compress/hold/expand evolution with the same original
Coulomb, radiation and alpha sources.

All 65 native fusion tests pass, including 13 new checks for scaled collision
geometry, separate adiabatic work, complete compress/hold/expand histories,
restart at a compressed volume, constant-history compatibility and atomic
domain/cancellation refusal. Independent 128-cell direct coupled Radau
integrations check temperatures, pressure, alpha moments and populations,
reaction accounts and all three work accounts, both with and without background
Coulomb exchange and radiation. Tightening the reference tolerance by 100 times
changes its energy accounts by less than `3.2e-9 J`; this checks temporal
coupling on that grid, not continuum energy-grid convergence. The runnable
256-cell example also completes, with its full energy residual below `6e-11 J`
in both cases. These checks build the actual native crates and dependencies;
the workspace DSR gate remains unavailable without DSR and sibling checkouts.

## Finite-time kinetic alpha heating in the light-ion network

`fusion::kinetic::advance_kinetic_fusion` couples the existing four-channel
thermal burn to a separately retained isotropic He4 energy distribution.
D-T and D-He3 births enter that distribution with their complete original
reaction-weighted product energy. Neither their energy nor their heat capacity
is assigned prematurely to the Maxwellian background. Canonical total amounts
still include their charge and baryons; thermal He4 stock is retained directly
through restarts rather than recovered by subtracting nearly equal inventories.
DD products and the D-He3 proton retain explicitly authored prompt deposition.
Alpha prompt-deposition fractions must be zero; the collision source determines
their electron/ion energy partition.

Each Maxwellian bath supplies the original finite-mass test-particle diffusion
law. A Scharfetter-Gummel birth/death generator satisfies separate discrete
detailed balance against `cell_width * sqrt(E) * exp(-E/kT_b)`. Positive
uniformization advances the frozen generator and both signed heat rewards with
the same bounded Poisson polynomial. Particles persist through their thermal
tail, with no absorbing zero-energy sink, terminal handoff, population clipping
or renormalization. Fixed-size trial storage supports 2–256 authored centers;
uniform-speed and logarithmic constructors preserve explicit endpoints. Birth
interpolation preserves number and mean energy. Both energy faces reflect.
Each channel's birth spectrum is collapsed to its original source event's mean
product energy; a resolved thermal/Doppler birth-energy spectrum is not supplied.

Midpoint coefficient refresh and symmetric reaction/collision composition use
step doubling over both thermal pools, all reaction accounts, each alpha bin,
alpha energy and both deposited-heat increments. This controls temporal error
on the authored grid, not energy-grid or physical-model error. Uniform speed
avoids the severe uniformization stiffness of tiny logarithmic low-energy bins.
An explicit upper-bin fraction screens inadequate high-energy support; it is
not a distribution-error bound. Grid refinement remains necessary, especially
for the width of an initially narrow birth distribution.

`KineticFusionState` publishes actual thermal energies, retained alpha energy,
total source amounts, and `p = 2*(Uh+Ue+Ualpha)/(3*V)`. It is not exposed as an
all-thermal gas state. Independent publication accounts check alpha number,
`delta(Ualpha)+Halpha_h+Halpha_e-Ealpha_born`, both original-source thermal
energies, every species, baryons, charge, and the full nuclear/escaping-energy
balance including alpha storage. Source objects, controls and distribution
are retained; cancellation or refusal publishes no partial result.

The original rigid-alpha checks include exact two-cell collision/reward
references and discrete Maxwellian refinement. Independent 128-cell Radau integrations check
seeded-alpha slowing and continuous DT births with finite thermal backreaction;
fixed-step refinement approaches that same-grid reference at second order.
Restart, all-channel Coulomb/radiation coupling, original-source and radiation
motion gates, positivity, cancellation and bounded collision-work retries are
also exercised. These native checks use the actual crates and dependencies;
the workspace DSR gate remains unavailable without DSR and sibling checkouts.

This rigid homogeneous owner requires explicitly dilute alpha number and energy
fractions, original He4/source identities, nonrelativistic grid support and
electron-Debye screening. Optional background Coulomb exchange excludes kinetic
alphas. Optional free-free radiation includes their ion charge only after a
separate fastest-represented-alpha velocity-variance check. No alpha self-collision,
anisotropic transport, confinement, fast-secondary reactions
or neutron transport is supplied. Original source/charge-reference admission
treats total He4 at the background temperature;
it does not certify classical occupation of an arbitrary authored alpha
distribution. Classical dilute test alphas are assumed, and the declared
number/energy limits do not quantify neglected alpha self-collisions. The frozen
thermal-reference leading-log approximation is distinct from a velocity-dependent
or sub-leading-order quantum
collision law. The thermal tail and separate two-temperature partition follow
the Fokker-Planck treatment discussed by Brown, Preston and Singleton,
<https://arxiv.org/abs/1106.1463>; their higher-order collision coefficients are
not claimed here. The `kinetic_fusion` example compares an explicitly prescribed
80% prompt electron split with finite-time alpha heating from the same sources.

## Prescribed fusion compression and separate thermal pressure work

`fusion::volume::advance_adiabatic_fusion_volume` moves either original D-T or
six-species light-ion plasma to an explicitly prescribed volume without
advancing time. Both monatomic thermal pools obey `p_s = 2 U_s/(3 V)` and
`T_s V^(2/3) = constant` at fixed populations. Mechanical work **on** each pool
is positive for compression and negative for expansion. Original source
energies independently check the separate work receipts; populations and source
identity are retained. Every possible source must cover the complete monotone
temperature path, including unpopulated products. Endpoint physical admission
suffices for this ideal adiabatic path: coupling scales as `V^(1/3)`, phase-space
density is constant, and thermal/rest-energy ratio scales as `V^(-2/3)`.

`FusionVolumeHistory` retains 1–256 original finite, positive volume knots at
strictly increasing times beginning at zero. Interpolation is explicitly linear
in log volume, giving constant fractional compression rate per segment. Holds
and authored endpoints retain their exact values, queries do not extrapolate,
and one knot describes a zero-duration history. Precomputed log ratios and
bounded binary search keep successful lookups allocation-free.

`fusion::network::advance_fusion_network_with_volume` couples that history to
the existing selected D-T, both D-D branches and D-He3 network. Every reaction,
Coulomb and radiation evaluation uses the actual stage volume and populations.
Exact adiabatic maps surround the midpoint source solve. Only the source solve
advances the external-heating clock; every accepted or attempted interval stops
at the next authored knot. The adaptive comparison includes both signed work
increments. An unrepresentable comparison midpoint refuses rather than
claiming zero local error. Smooth-segment refinement is second order; the
estimate remains numerical, not a global or experimental error certificate.

Moving states carry positive thermal energies independently of cumulative
transfer/work accounts. Compensated small increments retain tiny work; direct
positive scaling preserves the survivor after deep expansion. This avoids
recovering a small thermal energy as the difference of nearly equal large
initial-energy and work accounts. Publication checks both original caloric
energies, charge, baryons including emitted neutrons, and every species against
the reaction extents. The complete balance is
`delta(Uh+Ue) + En + Ephotons - sum(Qr*xir) - Eexternal - Wh_on - We_on = 0`.
Result objects retain the authored geometry, sources and numerical controls.
The rigid entry point remains available and reports exactly zero boundary work.

Fourteen new focused checks cover both source bindings, exact adiabatic scaling,
tiny volume changes, deep expansion, reversibility, source gaps, knot handling,
constant-history compatibility, atomic cancellation and physical-domain refusal.
The external-power benchmark solves `dU_s/dt = a U_s + P_s` independently and
checks fourfold error reduction on step halving. A separate direct-SI RK4
reference resolves the full compression/hold/expansion history, all five ion
amounts, four extents, eight product energies, both work accounts and thermal
transfers. Its 512/1024-step-per-segment refinement differs by less than
`5.5e-12` relatively in final temperatures and `4.7e-11` in extents; the net
electron work changes by `5.1e-7 J` after opposing work contributions cancel.
The runnable `compressed_fusion` example compares that history with a rigid
same-duration baseline using the same original source builder.

The ideal isotropic scaling is the fixed-ionization monatomic limit described
by Davidovits and Fisch, *Modeling turbulent energy behavior and sudden viscous
dissipation in compressing plasma turbulence*,
<https://www.osti.gov/servlets/purl/1432486>. Separate electron/ion pressure work
and reciprocal thermal exchange follow the two-temperature energy accounting
in *ALEGRA-HEDP v4.6*, SAND2004-5996, section 10,
<https://doi.org/10.2172/920747>. This owner prescribes homogeneous geometry;
it supplies no piston dynamics, shocks, anisotropic magnetic work, nonideal
correlation pressure or confinement. Charged products still thermalize
immediately, and selected neutron/photon escape remains explicit.

## Competing light-ion fusion and secondary-fuel burn

`fusion::network::advance_fusion_network` evolves explicitly selected D-T,
D-D→He3+n, D-D→T+p and D-He3 channels against one original six-species
`LightIonFusionPlasmaV1` bank. D-D products can feed secondary D-T and D-He3
burn during the same interval without trace seeds. Each independent D-D rate
receives the identical-pair event factor once and consumes two deuterons;
there is no additional assumed equal branching fraction. Every enabled
reaction must support the authored temperature interval, including initially
absent fuels. All four source domains intersect at 0.5–100 keV ion temperature.

Original shared thermal-event laws supply each channel's reactant kinetic
debit and both product energies. Per-product electron deposition fractions
are explicitly authored; all charged products remain and thermalize immediately,
and every neutron escapes. The source state supplies pressure, heat capacities
and sensible energy from actual populations. Nuclear baryons, charge and Q
replace chemical element/mass invariants for this owner only.

A fixed palindromic composition of exact binary finite-stock maps resolves
shared fuel competition at midpoint coefficients. Small withdrawals use
compensated populations; deep withdrawals use the analytic positive survivor,
without clipping or erasing representable traces. Shared Coulomb/radiation
kernels surround the reaction map, and held external powers enter its physical
clock once. Fixed-size trial storage allocates nothing on successful evolution.
Step doubling controls amounts, each channel extent, product energies, both
temperatures and thermal transfers. Smooth fixed-step refinement exhibits
second order; this is not an exact finite-step simultaneous network solution
or a global/physical error certificate.

Publication independently checks original source energies, every species
against channel stoichiometry, baryons including emitted neutrons, and charge.
The complete account is
`delta(Uh+Ue) + En + Ephotons - sum(Qr*xir) - Eexternal = 0`.
Results retain all sources, selected channels, deposition assumptions and
eight product-energy accounts. The original D-T API reuses the same thermal
kernels and reaction-energy law.

Twenty-four native gas-fusion tests pass, including the seventeen existing
D-T checks and seven new network checks. Independently converged RK4 references
cover initially pure D and mixed-fuel burn with Coulomb exchange, constant-Gaunt
radiation and unequal product-deposition fractions. Tests also cover exact DD
counting, fixed-step refinement, channel/source domains, conservation, volume
scaling, restart and atomic refusal. The compiled `light_ion_fusion` example
compares primary-only DD with secondary burn and demonstrates a mixed case
with the bounded Itoh radiation fit. The existing `dt_fusion` example retains
its numerical result. No confinement, fast-particle slowing, neutron transport,
opacity or relativistic product-spectrum model is inferred.

## Coupled fusion, Coulomb equilibration and escaping photons

`fusion::advance_dt_fusion_with_transport` adds explicitly selected
`DtFusionTransport` sources to the same finite-fuel, two-temperature burn.
Its original-bank Coulomb model refreshes actual fuel and ash populations;
its electron-ion bremsstrahlung model refreshes the charge-square density sum.
Selecting radiation explicitly assumes every emitted photon escapes this
optically thin parcel. The result retains both source models and reports
signed electron-to-heavy Coulomb transfer separately from outgoing photons.

Frozen reciprocal exchange integrates the temperature difference exponentially,
with reduced heat capacity `Ch*Ce/(Ch+Ce)`; it does not overshoot equilibrium
for a fixed coefficient. Freezing the instantaneous radiation coefficient
`A = P/sqrt(Te)` integrates `d(sqrt(Te))/dt = -A/(2Ce)` analytically.
Symmetric source composition surrounds the exact finite-fuel/external-power
step, and a predictor refreshes both coefficients at the interval midpoint,
including any selected temperature-dependent Gaunt factor. Capacities come
from each stage's actual composition. The method is second order for smooth
admitted evolution, as checked by fixed step refinement; the cutoff crossover
need not have a smooth local-error expansion. Step doubling remains an
estimate, not a global error certificate.

Driven reciprocal exchange also resolves its frozen relaxation exposure
`G*(1/Ch+1/Ce)*dt <= 1` before comparing trial paths. For
`D=Th-Te`, the frozen equation is
`D' = Ph/Ch-Pe/Ce-G*(1/Ch+1/Ce)*D`; saturated split maps can both return
`D=0` while missing the forced steady temperature difference, even when the
total energy balance closes. Differential imposed heating, selected particle
exchange, radiation, unequal confinement laws and nuclear/kinetic-alpha heat
activate this admission in the shared D-T/network source owner, including
kinetic and prescribed-volume callers. Zero conductance and unforced exchange
retain their unrestricted maps; equal imposed temperature rates alone do not
activate this restriction. Capacities are those at the start of the source map.
Changing populations still require coefficient refresh and temporal refinement;
this is not a bound on every intermediate relaxation time. A caller that
prohibits the required resolution receives the existing atomic `StepBudget`
refusal.

The publication account is
`delta(Uh+Ue) + En + Ephotons - Q - Eexternal = 0`.
Collisional exchange cancels internally, photon losses debit electrons once,
and neutron thermal energy remains included. Fixed-size trials and compensated
transfer sums retain allocation-free evolution and cancellation/refusal atomicity.
A restart reuses the retained operators on the actual new populations.

All seventeen focused native fusion tests pass, including seven coupled checks:
analytic no-fuel radiative cooling, reciprocal relaxation, two independently
converged constant-Gaunt RK4 burn references, second-order refinement, volume
scaling/restart, source-domain/cancellation refusal, and a separately converged
burn reference with the temperature/charge-dependent Itoh Gaunt fit.
Ten existing thermal-energy tests and thirteen ionization tests also pass.
The `dt_fusion` example runs both isolated and coupled cases and selects the
bounded Itoh fit with explicit cutoff/thermalization assumptions. Populated
D/T impose that fit's 157,900,000 K electron-temperature ceiling, independently
of the fusion reaction's ion-temperature support; no extrapolation is supplied.

## Finite D-T burn and two-temperature energy feedback

`fusion::advance_dt_fusion` evolves the existing
`FrozenIonizedGasStateV1` through the explicit L1
`DtFusionPlasmaV1` source binding. Its finite D+ and T+ inventories create
He4++ ash and an equal molar amount of escaped neutrons. Electrons retain their
amount. The original monatomic source bank evaluates the changed temperatures,
internal energy, pressure and heat capacities. Nuclear baryons and electric
charge replace chemical element/mass invariants for this owner only.

The declared closure is a homogeneous, rigid, quasineutral, weakly coupled,
nondegenerate, nonrelativistic Maxwellian parcel with separate ion and electron
temperatures. Every alpha is retained and thermalizes immediately, with a
caller-supplied electron/heavy deposition fraction; every neutron escapes.
Constant external powers are explicit. D-D and D-He3 reactivities exist in L1
but are not silently included in this D-T-only evolution.

Per-event ion removal is reaction weighted:
`E_rel = kT_i [3/2 + d ln(rate)/d ln(T_i)]`, plus `3kT_i/2` of
center-of-mass kinetic energy. Leading nonrelativistic two-body kinematics
partitions Q plus that energy between alpha and neutron. Thus escaping
neutrons carry their thermal contribution as well as their Q share. Nuclear
Q appears once, outside the source calorics' zero formation energy.
The final account independently checks source ion/electron energy, total
thermal change plus escaped neutron energy, baryons including neutrons, and
charge. These are numerical balances, not source-uncertainty certificates.

Each frozen-temperature reaction map solves the finite two-fuel depletion
analytically, including unequal fuel and surviving traces. An exponential
midpoint predictor with two-half-step error estimation resolves thermal
feedback. Fixed-size trials, compensated accumulated transfers, deterministic
math, bounded attempts and caller cancellation avoid per-trial allocation or
unbounded work. No negative population is clipped, trace seed inserted, or
equilibrium projection applied. Errors and cancellation publish no partial
state. Local tolerance does not certify global integration error.

The returned common state can enter another physical owner under an explicit
coupling schedule, subject to that owner's existing domain. The focused
regression suite includes independent RK4 burn endpoints, analytic depletion,
energy/charge/baryon/rest-mass accounts, volume scaling, temperature ownership,
zero-fuel and zero-time behavior, refusal/cancellation, repeatability/restart,
and an admitted handoff to the Coulomb thermal owner. That handoff uses cold
electrons to stay within the original classical Coulomb cutoff; it does not
extend classical transport to arbitrary hot fusion electrons.

Run the constructed source example with
`cargo run -p fs-gas --example dt_fusion`.
Its ideal monatomic source cards are explicit analytic inputs, not measured
caloric calibration. Spatial confinement and transport, MHD, opacity and
radiation transport, fast-particle slowing distributions, relativistic product
corrections, velocity-resolved quantum collision integrals, competing burn
networks, and neutron transport/deposition remain separate missing physics.

## Distributed mechanical work on native 1-D fields

With `experimental-euler`, `body_force::{BodyForce,ForceControl,ForceChange}`
serves both `dynamics::Euler1d` and `multicomponent::MixtureEuler1d`.
`apply_body_force` applies a held source without advancing the transport clock;
`step_with_body_force` composes exact source half-steps around the original
MUSCL/HLLC SSP-RK2 transport and its bounded retries. Uniform or per-cell
acceleration uses the density at EACH source stage; force density remains an
independent N/m³ input. One lazy flat scratch buffer is reused by each owner.
The all-zero source takes the original ordinary step without allocating it.

Represented momentum increment times midpoint velocity is the source work.
Species/mass and formation-energy references are unchanged by the source maps;
they add no heat. An impulse lost to momentum rounding cannot create energy.
`source_change` exposes compensated actual impulse, conjugate work and source
energy-representation residual. The complete step subtracts the work once and
retains the actual conservation defect without an energy correction. Existing
boundary enthalpy/pressure work, area force and atmosphere gravity stay distinct;
an additional body force must not reapply the bound gravitational acceleration.

Original EOS admission, explicit source-endpoint Mach limits and full-interval
CFL checks precede publication. A refused source, failed retry or cancellation
changes no accepted cell or clock. These are held-source stationary-mesh maps;
moving boundaries, finite networks, chemistry and molecular/wall exchanges keep
their existing separate composing owners. State/time-dependent source values
need temporal refinement; no new entropy-stability theorem is asserted.

### Original-source real-gas mechanical forcing

With `material,experimental-euler`, `RealGasEuler1d` also exposes
`apply_body_force` and `step_with_body_force` through the same represented
impulse/midpoint-work kernel. Every candidate retains the original
density-dependent NASA/Peng--Robinson calorics, pressure and sound speed;
neither an ideal-gas closure nor a constant heat capacity substitutes for it.
The selected native first-order/MUSCL SSP-RK2 transport and retries own the complete
force/transport/force transaction. Two lazy cell buffers retain the first
source endpoint and its original-source primitive cache, and are reused.
Source-endpoint Mach and full-interval CFL checks precede publication.

`step_with_walls_and_body_force` then advances the original finite-wall
heat/drag kernel over that accepted duration. This composition with walls
remains first order. Source impulse/work, wall impulse/heat/drive work and
boundary transport have distinct accounts; gas plus finite-wall storage
subtracts external body work exactly once. Cancellation or refusal publishes
neither cells, primitive caches, wall inventories, reports nor the clock.
Zero forcing preserves the original transport/wall path. The separate moving-wall
entry below composes this same source with prescribed boundary motion.
Hydrostatic well balancing, spatial viscosity and phase flash remain outside scope.
`tests/material_nonideal_flow/forcing.rs` checks dense-gas acceleration and
braking, gauge invariance, force-density/volume units, uniform forced transport,
endpoint-CFL retries, combined wall/body balances and late atomic cancellation.

### Distributed body work on moving intervals

`Euler1d::step_with_moving_walls_and_body_force` and its original mixture
counterpart compose the shared source maps with the existing extensive ALE
SSP-RK2 owner. The first half-source uses accepted volumes; the last uses the
admitted moved volumes and current densities. Uniform force density therefore
integrates the actual linear volume trajectory, while acceleration follows the
actual transported mass. Both source endpoints use relative-face CFL, original
calorics and the authored Mach bound.

Boundary pressure impulse and `p_star * wall_velocity` work remain in the
original end exchanges; independent body impulse, conjugate work and visible
energy-representation error remain in `source_change`. The complete balance
subtracts source work once. Each retry reconstructs both source halves and the
moving geometry from accepted fields, reusing the existing flat source buffer.
Cells, mesh, clock and accounting publish together after final admission.

Zero source preserves the original moving path, and zero wall velocity preserves
the stationary forced path. Constant area, affine mesh motion and reflecting ends
remain required; the frozen owner continues to refuse a bound gravitational
reference. Force values are held in original cell order. Chemistry, radial heat
transport and dynamically solved wall inertia retain their separate owners.

Six regression groups in `tests/body_force/moving.rs` define exact accelerating
translations with braking/reversal, original multicomponent calorics and absent
species, independent changing-volume impulse and boosted wall/body work, zero
identity, source/CFL retries and cancellation/refusal with unchanged physical
fields. Native execution remains outstanding in the constrained host session.

## Constant-pressure finite-rate chemistry

`chemistry::isobaric::IsobaricReactor` advances a sealed ideal-gas parcel at
prescribed pressure with total enthalpy conserved. It shares the original
extent integrator, backward Euler/SDIRK2 stages, adaptive comparisons, reduced
Newton factorization, source admission and conservation gates. The rigid U,V
entry points retain their original behavior. The prepared adapter reuses its
admitted storage across independent parcels; owning convenience results retain
the original mechanism and endpoint thermodynamic receipts.

Every trial computes its physical volume from nRT/p and evaluates kinetics at
the resulting actual concentrations. Reference-volume extent coordinates keep
stoichiometry exact; their rate Jacobian includes both volume and concentration
dilution derivatives. Temperature recovery uses original source enthalpy
increments, coefficient-combined reaction Delta h and Cp. Chemical formation
energy enters that single enthalpy constraint, never an additional heat source.
Entropy and elemental, charge and mass balances remain extensive at the
physical variable volume. Results expose initial/final V, H, U and outward
p Delta V, so the adiabatic first law can be checked independently.

The explicit-H entry point first resolves the initial thermodynamic state at
the authored composition. That resolved state defines initial volume and zero
boundary work; the temperature hint is not an additional energy source.
Unresolved energy targets, source domains, invalid pressure, cancellation or
bounded work exhaustion publish no accepted result. This is a homogeneous,
single-temperature ideal-gas reactor with prescribed constant pressure, not a
piston inertia, heat-transfer, nonideal or multiphase model.

Focused G1/G3 fixtures cover analytic isothermal dissociation, concentration-
sensitive bimolecular recombination, both stiff methods, full constrained
Jacobian differences, original reduced Newton work, actual binding-energy
changes and exact elemental energy-reference shifts. G4 exercises prepared
cancellation/reuse and zero-duration explicit-H semantics.


`fs-gas` owns L3 evolving gas systems. Constitutive and species laws remain in
`fs-thermochem`; this crate does not replace those sources with an independent
caloric model. Production dependencies remain pure Rust and Franken-owned.

## Material wall conduction in gas networks

With `experimental-euler,material`, `NetworkWallOperator::with_material_conduction`
binds original temperature or phase-dependent conductivity sources and explicit
planar length/area to each selected branch or plenum thermal contact.
`MaterialWallConduction::from_geometry` and
`from_phase_source_with_geometry` also admit the existing
`ReservoirLayerGeometry`: exact cylindrical or spherical radial resistance,
either contacted face, and declared angular area fractions. For uniform k,
the solid conductance is `2*pi*L*k/ln(r_outer/r_inner)` for a cylinder and
`4*pi*k/(1/r_inner - 1/r_outer)` for a sphere, scaled by the represented fraction.
The finite-wall trace radius represents its actual stored temperature; radial
area variation does not add storage or resolve transient through-wall diffusion.
`geometry()` returns the explicitly admitted multiplier and exposed face area;
original planar constructors retain their authored length/area and return `None`.
`length()` is the equivalent planar thickness `A_surface/C` for a radial binding.
This same binding feeds local walls, shared gas ports and wall-to-wall links.
The gas film and solid resistances add in series. Existing source-bound circular
gas transport still supplies the film from the actual gas mixture, density, T and p;
other contacts retain their originally authored film. The operator preserves
that separate film so repeated steps do not compound the same solid resistance.

In the default `WallConductionEvaluation::InitialState` mode, material k is
evaluated at the actual starting wall state of each source half, then frozen
for the original conservative implicit heat solve. The final
wall state must remain inside the original source domain. Each next half reads
the evolved wall temperature, including the existing reaction/wave composition.
Finite wall storage remains exclusively with the original contact and optional
material caloric binding; conductivity adds no second heat capacity or heat
supply. Any source refusal or late cancellation preserves all accepted gas,
wall temperatures and reports. Source lines, card identity and geometry remain
inspectable. The prepared temperature-only source refreshes without database
queries or allocation; phase-source refresh uses the original scalar resolver
with its selected claim pinned, retaining each query's evidence and conditions.

`MaterialWallConduction::with_evaluation(ImplicitEndpoint)` instead resolves
the original conductivity at every proposed caloric wall endpoint inside the
same bounded scalar heat equation. It requires the existing material heat-storage
binding and supports both temperature-dependent and phase-dependent sources.
The retained gas film remains separate and frozen over the source half. Candidate
storage, endpoint conductivity and actual represented heat must satisfy the same
conservative backward-Euler transaction. Increasing conductivity does not inherit
the fixed-conductance heat-flux bracket. Without an exact conductivity tangent the
solver uses bounded bisection; nonlinear conductivity carries no uniqueness claim.
Every source evaluation in the required heat bracket must be supported. Source,
work-budget and cancellation refusal preserve all accepted state.

`MaterialWallConduction::from_phase_source` uses `PhaseWallConductionSource`
and the same `MaterialWallCaloric::from_phase_curve` chart to derive BOTH T(h)
and liquid mass fraction. Thus equal temperatures on a melting plateau can
produce different source k(T,phi); no phase interpolation is invented by gas
code. The source card must match the chart, the storage chart must match the
conduction chart, and the original density must remain compatible with the
fixed reference mass/volume. Latent heat stays entirely in h. Missing phase
coordinates, a different chart, source extrapolation or an unavailable pinned
claim refuse before state publication. Temperature-only sources continue to
refuse energy-coordinate walls. Inspectors distinguish the two source forms.

`MaterialWallConduction::with_integrated_profile()` instead prepares the exact
original piecewise-affine conductivity line once and shares it across clones.
At every proposed gas/wall endpoint, the existing `fs-material`
Kirchhoff/film owner solves the material-face temperature and
`P = C * integral(T_face..T_wall, k(T) dT) = G_film * (T_face - T_gas)`.
The gas and wall retain their original caloric laws and exactly one heat debit.
The actual material profile must stay in source support; gas temperature may
lie outside the solid chart when the solved face remains supported. Exact
power derivatives fund Newton's slope, and the original source's maximum k
bounds gas-inversion error without inflating the relative physical heat budget.
The initial heat bound uses the integrated law's passive monotonicity. If an
unused storage endpoint has an unsupported material face, bounded source-aware
bracketing searches for a supported opposite-sign residual; it never clips a
face temperature or hides an unsupported physical root. The nested flux residual
is included in the final heat budget. Geometry preparation adds no per-query
allocation.

This opt-in supports local branch/plenum finite material contacts, prescribed
local reservoir contacts, and gas ports sharing one finite material inventory.
A prescribed reservoir retains its temperature and records the actual supplied
heat as external energy, with zero wall storage. Shared ports solve their gas
endpoint equations against one trial wall state and debit the retained wall
only by the compensated sum of accepted heats. The shared Newton derivative
uses the original power partials and each gas's original heat capacity, never
the nonlinear profile's secant conductance. Every source query remains in the
existing shared work budget; authored gas films are never overwritten.

At a trial shared-wall temperature, the artificial pair with the initial gas
temperature may have an unsupported face despite a supported physical root.
The inner gas solve starts at the original gas-domain temperature nearest the
wall, then searches for a supported opposite-sign endpoint. Both inner gas
and outer wall searches retain original source refusals and bounded work.
Later wall-link stages re-admit actual final gas-facing material profiles
without applying a second heat transfer. Wall links themselves still refuse
this opt-in until a two-layer endpoint law is supplied. Phase-dependent
k(T,phi) also refuses: one lumped wall h cannot provide a spatial phase field. Selecting
`InitialState` after opting in refuses at binding rather than silently replacing
the integrated law with a point sample. Source-domain or cancellation refusal
preserves the accepted network and all original wall inventories.

The initial-state mode retains constitutive lag. Ordinary implicit endpoint
evaluation closes that local lag but still samples one wall temperature. The
integrated profile adds a quasi-steady through-wall constitutive solve; all
modes retain fixed geometry, lumped-wall storage, backward-Euler time integration
and frozen gas-film/drag coefficients. None resolves transient internal diffusion.
For centre-temperature storage, the conduction length is the surface-to-centre
distance. Adiabatic contacts refuse. Original phase-aware spatial and mechanical
owners remain responsible for moving fronts, wall support, shrinkage and liquid
motion; a thermal wall cannot remain a mechanical support after melting merely
because this heat-transfer operator still exists. Source endpoint admission is
not a continuous source-domain certificate between sampled states. Focused
G1/G3 tests check independent two-body heat balances, source-only conductivity
and capacity swaps, combined gas-film/solid response, actual chemical-rate and
pressure feedback, endpoint-domain refusal, cancellation and deterministic replay.
Phase regressions additionally check equal-T distinct latent inventories,
sensible-to-latent crossing, three-rung time refinement against an independently
integrated latent-plateau energy equation, k refresh without compounding resistance, exact
original-chart matching, pinned-source refusal and native reacting-network
energy/rate feedback. These manufactured fixtures are not measured-material
or wall-geometry qualification.

## Experimental frozen-composition Euler flow

Feature `experimental-euler` exposes `dynamics::{FrozenGas, Euler1d}`. It is an
Estimated engineering baseline for one-dimensional, inviscid, stationary variable-area flow
of a single frozen ideal-gas mixture. It advances actual compressible waves,
contacts, shocks and rarefactions without the low-Mach lattice assumption.
It does not grant the separate entropy/invariant-domain lane's theorem claim.
This implements the first transport slice of
`frankensim-ext-gas-1d-network-vvnx`; that broad network bead remains open.

### Source and thermodynamic state

The immutable `FrozenGas` owns the complete original NASA-9 mixture and a
caller-declared temperature interval. All active species must support that
interval without gaps. Crossing source regions requires the original source's
explicit caloric-continuity construction. Bounded directed-rounding interval
subdivision of each component's original Cv/R polynomial admits positive Cv
over the whole path. Failure to establish that bound refuses construction.
The subdivision uses a fixed-size stack and at most 4095 work items per region.

Pressure is `rho R T`; frozen sound speed is `sqrt((Cp/Cv) R T)` at the same
temperature and composition. Temperature inversion uses the original source's
reference-invariant internal-energy difference and analytic Cv, with safeguarded
Newton/bisection and explicit absolute/relative energy tolerances and iteration
limits. There is no constant-gamma substitution for a variable-Cp source.
Full source-bearing evaluations remain available through `thermodynamics`.

Conserved cells contain mass density [kg/m³], axial momentum density
[kg/(m² s)] and referenced total-energy density [J/m³]. The energy variable is

`rho * (e(T) - e(Tmin) + u²/2)`.

The original signed `e(Tmin)` remains in the owner. Physical total energy is
recovered by adding `e(Tmin)*rho`, and physical energy flux by adding the same
reference times mass flux. This exact linear change of conservative variables
avoids rounding away sensible energy under large formation references. It is
valid for the declared frozen composition and cannot silently transfer to a
reacting or composition-changing cell. No positivity requirement is imposed on
the original absolute formation/internal energy.

### Conservative transport and boundaries

The numerical flux is HLLC using both states' actual frozen acoustic speeds.
Contact intermediate states must have positive density/pressure and internal
energy inside the admitted source interval. Degenerate or inadmissible contact
states use an HLL fallback; its intermediate state is checked, with bounded
wave-speed expansion if required. These local checks are not a proof of global
entropy stability or positivity for arbitrary EOS data.

`NormalConserved` and `FrozenGas::hllc_normal_flux` expose the same closure to
two-dimensional spatial owners in an orthonormal normal/tangent frame. They
carry both momentum components and full referenced total energy. Transverse
momentum and kinetic energy follow the same HLLC/HLL waves as the normal
equations; HLL intermediate-state admission includes heat produced by numerical
shear mixing. No extra constant-gamma source or independent passive-scalar flux
is introduced. `normal_conserved` refuses loss of sensible-energy precision
when adding transverse kinetic energy, and `normal_primitive` removes that
energy exactly once before the original caloric inversion. The stationary slip
wall adapter reuses the thermally perfect shock/rarefaction wall solve and
returns zero mass, transverse momentum and energy flux. Spatial reconstruction,
mesh evolution and timestep acceptance remain with the consuming spatial owner.
Focused G0/G3/G4 tests cover axial correspondence, moving shear contacts,
transverse-boost covariance on both flux branches, shear-heating source-domain
refusal, slip-wall pressure and cancellation. This adds no entropy theorem,
viscous wall model or physical shock-thickness claim.

`NormalPrimitive` and the `*_normal_flux_from_primitive` entries admit
reconstructed density, normal/tangential velocity and pressure through the
same source and conservative-energy precision checks. They reuse the original
source-computed thermodynamics in the unchanged HLLC/HLL and wall kernels,
avoiding an immediate inverse temperature solve on energy just constructed
from those same primitives. Callers cannot supply a cached sound speed or
substitute a caloric model. Two added regressions compare conservative and
primitive entries for contacts, supersonic flow, HLL shear mixing and both wall
branches, plus source/precision/cancellation refusals. All 154 focused native
gas library tests pass on the pinned nightly.

`SpatialConserved`, `SpatialPrimitive` and `SpatialFaceFlux` extend that same
closure to a normal and two orthogonal tangent directions. The state carries
`rho * (e(T)-e(Tmin) + (u_n²+u_t1²+u_t2²)/2)`. All three momenta and the full
energy participate in one HLLC/HLL solve. A private const-generic kernel is
shared by the planar and volumetric entry points, including star-state
admission, bounded HLL wave-speed expansion, and physical flux evaluation;
it does not fork the Riemann model or append independent scalar advection.
The three-dimensional constructor, decoder and primitive fast path preserve
the same original source and caloric-precision admission. Wall pressure uses
the existing thermally perfect shock/rarefaction closure, with zero mass,
both tangential fluxes and work at a stationary adiabatic slip wall.

The spatial face tests check bitwise planar reduction when the second tangent
vanishes; double-shear contacts; rotations and Galilean boosts of both tangent
components on HLLC and HLL branches; combined shear-heating domain refusal;
tangentially invariant wall pressure; source-owned primitive/conservative
agreement; and precision/domain/cancellation refusal. This is the closure for
the multi-D owner in `fs-flux`, not a three-dimensional mesh integrator here.
The Estimated claim and the entropy, vacuum, wall-model and physical-validation
limitations remain unchanged. All 162 focused native gas library tests passed
on the pinned nightly through the original dependency sources; the complete
workspace DSR gate was unavailable and is not implied by this narrower result.

`moving_hllc_normal_flux`, `moving_hllc_spatial_flux` and their reconstructed
primitive entries expose conservative transport through translating internal
faces. They reuse the original HLLC/HLL problem in the face's normal rest frame
and return laboratory mass, all momenta and full referenced energy flux. The
existing one-dimensional moving-volume solver now uses this same transformation.
The source reconstructs face-frame caloric storage; no subtraction of large
laboratory kinetic energies or independent shear advection is used. Its HLL
middle-state check includes both tangential kinetic energies. The transformation
uses fixed-size stack arrays; source evaluations retain their original limits.

The returned flux is `F - w_n U`, and its signal bound is relative to the face.
Laboratory momentum gains `w_n * mass_flux` on transformation from the face
frame; laboratory energy gains `w_n * normal_momentum_flux + w_n^2/2 * mass_flux`.
A comoving pressure/shear contact therefore transfers zero mass and tangential
momentum while retaining pressure traction and `p*w_n` work. Zero normal speed
returns the original stationary result exactly. These are the moving-face
relations in Springel (2010), [arXiv:0901.4107](https://arxiv.org/abs/0901.4107).
The consuming spatial owner still supplies consistent moving volumes, opposite
shared-face signs, geometric conservation and time integration. No general
moving-mesh solver, entropy theorem or new material qualification is implied.
G0/G1/G3/G4 face regressions cover stationary/axial/planar reduction, Reynolds
transport, double-shear contacts, the comoving slip-wall limit, full momentum
and work covariance on HLLC/HLL, source-frame precision and cancellation.

MUSCL uses monotonized-central primitive slopes, reduced together when needed
to retain positive density/pressure and the declared source interval. SSP RK2
uses one shared oriented flux per face. Periodic end faces share the same flux
bits; reflecting stationary walls impose exact zero mass and energy flux;
transmissive boundaries use zero-gradient ghost states and are not claimed to
be nonreflecting subsonic outlets.

Every stage checks its actual face-speed CFL against the caller's value, at most
0.45. Failed trials halve dt within an explicit retry cap. CFL-only retries
retain the requested spatial reconstruction; actual admissibility failures use
a first-order retry. Every accepted stage and final cell must remain in the
source domain. No conservative variable, temperature or pressure is clipped.
`StepReport` retains the actual time increment, maximum stage CFL, fallback
counts and extensive mass/momentum/energy balances against the two RK stages'
integrated boundary fluxes.

Workspace is allocated once, with a hard cap of 1,048,576 cells and fallible
allocation. Successful steps allocate no numerical workspace. Cancellation is
polled through cells, faces, source admission, inversion and retries. A refusal
or cancellation preserves the complete original cell field and time. Operations
are ordered deterministically on one arithmetic target; cross-ISA bit identity
and parallel scaling have not been measured.

### Focused native evidence

Eight tests of the complete source compiled and executed through the actual
pure thermochemistry dependency closure. They cover source-gauge invariance
under ±1e18 formation offsets, variable-Cv/high-speed EOS round trips, HLLC
stationary contacts, uniform state and reflecting-wall conservation, strong
rarefaction positivity, source-domain refusal and early/mid/late cancellation.

Measured regression results on the authored fixtures:

| Case | Resolution | Result |
| --- | --- | --- |
| Smooth contact advection, density L1 | 40 / 80 / 160 cells | 7.98344e-4 / 1.96127e-4 / 4.72099e-5 |
| Sod tube, independent exact gamma-law oracle, density L1 | 120 / 240 cells | 4.01893e-3 / 2.20519e-3 |
| Lax tube, same independent oracle, density L1 | 120 / 240 cells | 1.48205e-2 / 9.28277e-3 |
| Variable-Cp acoustic wave | 128 cells | Measured speed 1.15993219 vs source speed 1.15966702 m/s |

The independent oracle uses shock Hugoniot and rarefaction-invariant relations,
not production HLLC/source-inversion code. A nonlinear variable-Cp shock also
executes with conserved boundary accounting. These are manufactured/conformance
tests, not experimental validation or a performance benchmark.

`cargo run -p fs-gas --features experimental-euler --example shock_tube -- 240`
runs a one-metre N2 tube with the 200–1000 K NASA-9 coefficient region already
retained in the thermochemistry fixtures, and emits the resulting field as CSV.
The native run advanced 292 steps to 4.739365351855e-4 s, with zero reported mass
balance residual and -7.68e-15 J referenced-energy balance residual. The material
association is explicitly declared by the example; this is not a validated
nitrogen shock-tube experiment.

### Scope boundary

No evolving species transport, chemistry/flow splitting,
junctions, valves, viscous stresses, heat conduction, general moving geometry, vacuum state,
real-gas EOS, multidimensional flow, adjoint, entropy certificate or experimental
validation is supplied by this Euler slice. A source-domain refusal is an actual
refusal and does not extend a NASA fit to hotter or colder shocked gas.

References: Toro, Spruce and Speares, *Shock Waves* 4 (1994), 25–34,
DOI 10.1007/BF01414629; Ketcheson, LeVeque and del Razo, *Riemann Problems and
Jupyter Solutions* (2020), Euler and approximate-Euler chapters.

## Stationary area and nonuniform-mesh extension

`Mesh1d` retains ordered physical face positions, positive face areas, and
integrated cell volumes. `Euler1d::new_on_mesh` consumes that geometry without
resampling or substituting a midpoint volume. The uniform constructor remains
available. `cell_width()` and `cross_section()` return `Option` so a nonuniform
or tapered field cannot advertise a nonexistent uniform dimension.

MUSCL slopes use the actual distances between cell centres. Each physical face
flux is multiplied by its area exactly once and shared between neighbouring
volumes. The pressure-area momentum source is paired algebraically with the
pressure traction, preserving a resting uniform-pressure state even across
density contacts in a taper. Stationary area changes do no energy work.

Both stages enforce a local area/volume-scaled face CFL, including any HLL
wave-speed expansion. `StepReport::area_momentum_impulse` records the geometric
force integral. The residual subtracts both boundary transfer and this impulse.
Periodic boundaries require equal end areas; they share the same flux bits.

Conservative state construction also checks that stored momentum and energy
retain the requested sensible energy within the existing caloric tolerance.
Finite values alone do not admit a state whose kinetic energy has rounded away
its thermal energy.

Fifteen focused native gas tests pass, including the original eight tests and:
nonuniform-taper rest with density contacts; finite-amplitude tapered pulse
mass/energy/momentum closure; mesh and periodic-boundary refusals; cancellation;
caloric-round-trip refusal; and a moving variable-Cp Mach-3 shock whose field and
reports are bitwise invariant under a +1e18 formation-reference shift.
The latter checks front position against independent Rankine-Hugoniot motion
and checks shock heating and entropy increase from the original source.

For smooth periodic contact advection on a stretched mesh, density L1 error
falls from 1.66346e-4 at 40 cells to 3.67860e-5 at 80 cells (4.52-fold).
This is a manufactured convergence result, not a throughput benchmark.
Source geometry/quadrature accuracy and the quasi-one-dimensional smooth-area
model remain explicit assumptions; abrupt junction losses are not supplied.

## Gravity-balanced isothermal atmospheres

`Euler1d::new_isothermal_atmosphere` adds signed constant axial gravity to the
existing frozen-mixture Euler owner. Its `IsothermalAtmosphere` contains typed
acceleration, temperature, pressure and position. The source's original gas
constant determines the exponential hydrostatic profile; original NASA calorics
determine energy and sound speed. Nonuniform constant-area cells receive
analytically integrated density/energy averages, with a regular small-gravity
series. Unsupported temperatures, unrepresentable profiles, nonconstant area,
inconsistent volumes and periodic boundaries refuse before publication.

Both first-order and MUSCL reconstruction operate on perturbations about the
retained reference. Source-evaluated reference faces are prepared once. Pressure
traction and the reference gravitational force cancel before the update; the
remaining density perturbation receives its actual gravitational force. The
gas evolves adiabatically after initialization: the reference is not a
thermostat and does not freeze disturbances or substitute constant gamma.
Successful steps reuse the original RK workspace and allocate no new numerical
arrays. Original source-domain admission, CFL checks, retries and cancellation
remain in force.

The shared one-, two- and three-dimensional HLLC kernel returns the physical
flux directly for identical admitted traces and the exact pressure-only flux
for a stationary pressure contact. This avoids roundoff transport from a
redundant star-state multiplication/division, including across a stationary
tangential shear jump. Ordinary shocks, moving contacts and HLL fallback still
use the original source-admitted wave solver.

Gravity's energy work uses the same face mass flux as continuity, making the
discrete sum of gas energy and `rho*phi(cell_center)` conservative up to boundary
transport and rounding. `StepReport` exposes the gravitational impulse and
work separately from area force and includes both in the balance residual.
`gravitational_potential_energy` exposes the retained midpoint potential
quadrature; this is not an exact within-cell gravitational-energy integral.
The focused G1/G3/G4 cases cover exact atmospheric averages, stationary balance,
disturbance evolution with closed energy accounting, analytic free fall under
refinement, formation-reference/coordinate covariance and cancellation rollback.

This reference-state method is bounded to a frozen ideal mixture, constant
acceleration and constant area. It does not claim balance for arbitrary
nonisothermal atmospheres, self-gravity, multicomponent/reaction transport or
an entropy/positivity theorem. Large disturbances may still refuse source
admission. Frozen networks use the explicit shared-potential constructor described
below. Arbitrary caller-owned end ports still refuse without their own declared
gravitational-potential and energy convention.

## Isentropic hydrostatic atmospheres

`Euler1d::new_isentropic_atmosphere` extends the same gravity-balanced owner to
uniform-entropy references with nonuniform temperature. It solves
`h(T)-h(T_ref) = g*(x-x_ref)` using the original positive-Cp NASA source and
recovers pressure from its source entropy difference. The lapse is `g/Cp(T)`;
neither a constant gamma nor a prescribed linear temperature fit replaces the
source. Internal-energy and entropy differences cancel their original reference
constants before arithmetic. `hydrostatic_reference()` distinguishes both
families; the existing `atmosphere()` accessor remains specific to isothermal
references.

Source-evaluated face states are retained once. Conserved density and referenced
energy are integrated with a cancellable composite four-point Gauss rule, using
`AtmosphereQuadrature`'s relative convergence tolerance and explicit panel cap.
Successive-rule differences are estimates, not certified integration bounds.
Caloric inversion retains the source owner's work cap and tolerance, tightening
the reference solve toward floating-point accuracy where the requested tolerance
allows. Source-domain escape or exhausted quadrature/inversion budgets refuses
before returning an owner. This avoids confusing spatially averaged primitive
fields with the decoded state of conserved cell averages.

The same prepared reference, perturbation reconstruction and shared-flux gravity
work preserve the discrete atmosphere and evolve real disturbances. There is
no added thermodynamic solve, integration or allocation in time stepping. Physical
gravity-force consistency converges with the reference quadrature and caloric
solve. G1/G3/G4 cases check analytic constant-Cp polytropes, independent 70-digit
variable-Cp conserved integrals, force balance, gas-plus-potential energy, formation
reference invariance and cancellation/restart. The existing constant-area,
constant-acceleration, frozen ideal-mixture and physical-boundary restrictions
apply; network coupling requires the shared-potential declaration below. This is a convectively neutral source isentrope, not arbitrary
entropy stratification, an evolving chemistry/phase equilibrium, or a global
entropy/positivity theorem.

## Physical reflecting-wall Riemann boundary

Reflecting ends now solve a one-sided thermally perfect Riemann problem
directly. Incoming gas connects to a stopped wall state by the source-owned
shock Hugoniot; the velocity jump is factored from the root equation so weak
shocks do not converge to an unchanged spurious root. Retreating gas connects
to the wall by an isentropic rarefaction, integrating Cp(T)/c(T) in sqrt(T)
coordinates and recovering pressure from the source entropy increment.
Quadrature is bounded, cancellable, and split at original caloric seams.

The wall imposes exactly zero mass and energy transfer and positive normal
pressure. It does not reuse an HLL ghost-state momentum flux that can become
artificially tensile during expansion. Unsupported hot or cold wall states
refuse; the source is not extended to an unadmitted vacuum. The actual wall
wave and acoustic speeds enter the same stage CFL calculation.

Twenty-one focused native gas tests pass. Six added checks compare weak and
strong constant-Cp shock/rarefaction pressures with analytic relations, check
variable-Cp Hugoniot and entropy balance, compare expansion with an independent
2,000-panel quadrature, exercise orientation/formation-reference invariance,
and verify source-window refusal and cancellation. These qualify the
one-dimensional inviscid boundary implementation, not viscous wall layers,
heat transfer, physical shock thickness or a fully discrete entropy theorem.

## Adaptive sealed finite-rate reactors

`chemistry::{advance_adiabatic, advance_adiabatic_energy}` advances an actual
caller-declared duration at fixed volume and total internal energy. The former
obtains the energy from the original feed temperature; the latter accepts the
actual stored energy from another conservative owner. The complete original
`IdealGasKineticsV1` mechanism, source receipts, canonical species amounts,
elapsed time, final pressure/temperature and balances remain in the result.
This implements a bounded homogeneous reactor slice of
`frankensim-ext-gas-chemistry-ladder-paqh`, not the complete combustion ladder.

The solver advances reaction extents with damped backward Euler (the default)
or optional second-order L-stable SDIRK2, both with adaptive step doubling.
Each candidate uses the last accepted composition plus local reaction extents
through the exact admitted stoichiometric matrix. The origin is refreshed only
after acceptance; the original feed and fixed total energy remain untouched for
whole-call mass, elemental, charge, energy and entropy audits. No species clipping,
positive trace seeding or equilibrium projection is introduced. Temperature solves
the original total internal energy at each candidate composition, with the
formation contribution present exactly once. No separate reaction heat source
is added. The source-path/positive-Cv interval admission is the same L1 helper
used by the Euler EOS, including possible products absent from the feed.

Caloric roots use `sum(c_i * Delta u_i(T_initial,T))` plus the compensated
original reaction energy accumulated through accepted extents, instead of
subtracting two gross formation-energy totals. Source reaction energies are
combined before polynomial evaluation by the same L1 kinetic owner. This
preserves actual reaction binding energies and variable Cv while canceling
conserved reference offsets. Rejected stages do not change the accumulated
chemical energy. The prepared workspace retains one species energy vector and
reaction energy/entropy vectors; nonlinear trials allocate no additional storage.

The physical `energy_residual_j` uses these source-integrated differences.
`absolute_energy_residual_j` independently subtracts the original endpoint
source energy from its target. Each accepted species update also records the
original-energy equivalent of `c_new-c_old-nu*xi` as
`composition_reference_roundoff_j`. It is a numerical inventory diagnostic,
not heat. The independent absolute-energy check allows this measured term and
a magnitude-scaled floating-point allowance; this is not an interval certificate.
Element, charge and mass gates retain their original tolerances. Explicit
absolute energy targets whose floating-point spacing exceeds the requested
Joule tolerance refuse; temperature-initialized evolution retains its precise
source-relative equation. G1/G3 tests independently recover temperature from
the published product amounts, including changed binding energy and molecular
count under positive and negative common elemental references.

The Newton Jacobian combines the analytic concentration and temperature rate
derivatives with `dT/dxi_j = -sum_i(nu_ij u_i)/sum_i(c_i Cv_i)`. A bounded scaled
partial-pivot solve and positivity/energy backtracking support stiff reactions.
Stage, Newton and root scratch is allocated once and reused. Local error control
compares the complete step with two half steps, covering species and temperature;
the two-half-step result is accepted without extrapolating species. Backward
Euler is first order; SDIRK2 is second order for smooth homogeneous reactor
ODEs. Both use estimated local errors rather than global accuracy certificates.

Caller controls separately bound absolute energy, absolute amount, relative
species/temperature error, relative element/charge/mass balance, numerical entropy
decrease, minimum/initial/maximum time step and all nonlinear/root/attempt work.
Every accepted step, SDIRK stage and complete integration check actual
ideal-mixture entropy using compensated source-relative increments. Original
reaction entropy coefficients cancel conserved elemental references before
evaluation; genuine species entropy differences still change reaction affinity.
Analytic Cv/T integrals, stable changes in c ln(c), the exact zero-population
limit and actual volume expansion supply thermal and mixing contributions.
Fixed-pressure evolution therefore retains both changing molecular count and
expansion entropy. No large absolute source entropies are subtracted for the
second-law gate. Independent isomer-mixing and dissociation formulas cover both
time methods under positive and negative elemental entropy-reference shifts.
Heating and cooling checks independently include the reservoir entropy change;
rebasing accepted extents preserves the same parcel entropy and temperature.
Final element, charge, molecular-mass and energy balances are checked again
before publication. A budget exhaustion, unsupported source state, arithmetic
failure or cancellation returns no partial integration. Zero duration preserves
the admitted initial state. No ignition delay or empirical rate accuracy is
inferred from successful integration of a caller-authored mechanism.

The focused reactor fixtures cover analytic reversible-isomer time evolution,
stiff relaxation to the independent Gibbs/UV endpoint, exothermic energy and
entropy balance, molecular-count-changing dissociation, extensive scaling,
one-correction Newton convergence, overflow/refusal and cancellation. Source
interval admission also rejects negative interior Cv despite positive endpoint
values. These are numerical/constitutive conformance tests, not experimental
validation of a real-fuel mechanism.

### Ionized species and numerical charge conservation

The same homogeneous reactor accepts charged mechanisms admitted by the
original `fs_thermochem::kinetics::IdealGasKineticsV1` owner. Integer stoichiometry must
conserve both original elemental and charge vectors, and the declared source
molar masses must agree with every reaction. Ions and electrons contribute to
the original species pressure and internal energy. Ionization energy belongs
to those source formation references and is never added as a second heat term.
For reversible channels, the original source Gibbs energies determine reverse
rates and the concentration-normalized equilibrium constants.

Exact stoichiometric conservation does not by itself bound floating-point
charge drift. Before each accepted step is rebased, the reactor compares its
actual concentrations against the untouched initial concentrations. Final
publication also compares the actual returned mole amounts against the input
mole amounts, including the effects of volume conversion. Both checks use
`IsochoricKineticsSolve::relative_conservation_tolerance` and the residual

`abs(sum_i z_i (n_i - n_i_initial)) / max(sum_i abs(z_i) n_i, sum_i abs(z_i) n_i_initial)`.

An identically zero carrier inventory has zero residual. Charges and carrier
amounts are scaled before multiplication, and compensated sums preserve the
canonical species order. Neutral species never enter the denominator. The
audit takes linear work in species count and allocates no storage.
`IsochoricKineticsBalance::maximum_relative_charge_residual` retains the largest
accepted-step or final-publication residual. An excessive or nonfinite residual
refuses the integration with `IsochoricKineticsError::Balance`; the solver does
not clip, renormalize, neutralize or publish an inconsistent endpoint.

The conserved target is the authored feed charge, not necessarily zero. This
is a homogeneous constituent/energy evolution law: admitting a nonzero charge
does not establish an electrostatic equilibrium or supply its field energy.
All species share one source temperature and ideal-gas convention. Separate
electron and heavy-particle temperatures, electric-field acceleration, ambipolar
diffusion, non-Maxwellian electron distributions, radiation-driven ionization,
and nonideal plasma corrections require separate owners and are not supplied
by this reactor.

Both owned `IsochoricKineticsResult` and borrowed `IsochoricKineticsView` expose
`frozen_state`. This explicit owning conversion retains the accepted actual
amounts, source bank, volume and common temperature in L1's
`FrozenIonizedGasStateV1`, for nonequilibrium charged-transport consumers. It
performs no Saha/equilibrium projection. The original reactor still preserves
an authored nonzero feed charge, but conversion to this quasineutral material
lane refuses it. State cloning shares the immutable sources and populations;
the prepared reactor's ordinary advance remains allocation-free.

The ionization regressions use `H + e- <-> H+ + 2 e-`, ideal monatomic
Sackur--Tetrode thermodynamics and NIST's H I ionization energy. They compare
detailed balance and the final adiabatic state with independent Saha formulas,
and finite-time evolution with an independent scalar integration. They check
electron pressure, ionization cooling, mass/energy/charge balance, preserved
nonzero feed charge, and refusal of trace-carrier drift over extensive scales.
Exactly zero initial electrons remain zero for this electron-impact channel;
no artificial seed is inserted. Its forward rate is explicitly synthetic, so
these numerical tests do not validate an empirical ionization time or discharge
mechanism. Source reference: NASA RP-1311, Part I, section 3.7,
https://ntrs.nasa.gov/citations/19950013764, and NIST H I ionization energy,
https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=H+I.

An additional source-backed regression runs the original stiff SDIRK2 reactor
with L1's NIST H BEB electron-impact law, including its Maxwellian temperature
derivative and the original Gibbs-derived reverse law. Ionization increases
the electron count and consumes thermal energy while mass, nuclei, charge and
total internal energy remain conserved. This proves that the retained physical
rate reaches the existing evolution path; it is not empirical validation of
ionization time or an electrical discharge.

## Frozen electron/heavy energy evolution

`chemistry::energy::advance_frozen_energy` evolves the actual source-bound
`FrozenIonizedGasStateV1` at fixed species amounts and volume, with independent
constant electron and heavy-particle powers over the requested interval.
Positive power deposits energy and negative power removes it. An absent
electron population cannot accept electron power or acquire a seed. Electrical
consumers supply their actual, separately resolved carrier work; this owner
does not invent a field or add total Joule work a second time.

Elastic energy exchange uses the same source-retained electron-neutral
momentum cross sections as material drag:
`G = 3 R n_e sum_target[N_target (m_e/m_target) k_m(Te)]`, with electron amount
`n_e` in mol, neutral number density `N_target` in m^-3 and
`Q_e_to_h = G (Te-Th)`. Every present neutral requires exactly one original
momentum source. The declared target/electron speed limit also bounds
`sqrt(m_e/m_target)`, so cold but unrealistically light targets cannot pass a
massive-target approximation. Charged fraction, source temperature intervals,
cross-section coverage and numerical work remain explicit. These limits do
not establish that neglected Coulomb or inelastic processes are insignificant
for an arbitrary experiment. The model is Banks (1966), Eqs. 16–20,
https://ntrs.nasa.gov/api/citations/19660020172/downloads/19660020172.pdf.

Each first-order substep freezes G at its starting state and solves a monotone
implicit scalar exchange. The same signed exchange enters the two pool energy
targets with opposite signs. Full/two-half comparison controls an estimated
local error, including the changing collision rate; it is not a global error
certificate. Every caloric inversion uses original-source internal-energy
increments relative to the original temperatures, canceling formation offsets
before rounding. Energy and temperature residuals both gate convergence, so
a trace electron pool cannot hide a large Te error under a parcel-wide energy
tolerance. Nonlinear trials allocate no heap storage. Source gaps, invalid
domains, exhausted work and cancellation publish no partial endpoint.

The result retains the initial and final source-bound states, original
collision sources, original controls, separately deposited energies, signed
exchange and actual per-pool energy residuals. This is a frozen-composition
thermal step. It does not supply two-temperature molecular detailed balance,
ionization population evolution, photon absorption, field energy or charged
spatial transport; those owners must exchange their own conservative sources.

## Common material gas integration

Optional feature `material` exposes
`material::advance_closed_parcel_adiabatic`. It evolves an existing
`fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel`, using its
actual stored internal energy as the one target. The complete product bank
must match the mechanism's original NASA-9 models and retain every active
feed transport card, source-pack identity and declared water identity.
This adapter refuses mechanisms containing charged species before integrating:
its Wilke viscosity and Wassiljewa--Mason--Saxena thermal-conductivity closure
belongs to neutral dilute gas. Source identity alone does not make that closure
an electron/ion transport model. A charged-mechanism refusal leaves the
original parcel unchanged, including for a zero-duration request.

Generated amounts pass through `ThermochemicalGasModel::with_mole_amounts`,
which uses the existing exact-simplex closure and refuses unrepresentable
positive traces. The existing material composition/energy replacement resolves
all final transport and acoustic properties at the preserved mass and volume.
The final material state is checked against both the original stored energy
and the kinetic endpoint, including mass, density, temperature, pressure and
composition-roundoff discrepancies. Mass and volume retain their original
bits; zero duration retains the entire original parcel. Publication follows
the last source, balance and cancellation check.

Five material-consumer regressions exercise finite elapsed reaction against an
analytic solution, reaction-driven heating and pressure rise, updated common
transport/acoustics, exact storage preservation, source/water identity refusal,
zero-time identity, early/late cancellation, exhausted work and final transport
domain refusal. The numerical reactor and source model are the production
implementations; the bounded native material test closure selects original
material storage/resolver/caloric bodies and excludes unrelated runtime and
mechanical modules. It does not establish a full `fs-material` or workspace
build. The existing workspace/runtime resource limits remain applicable.

Across the final focused native runs, all eight reactor tests, fifteen Euler
tests and five material-reactor tests passed. Eighteen selected original
material tests also passed, including the two generated-amount regressions;
the complete original thermochemistry library passed all 111 tests. The
material-consumer run compiled the actual `fs-gas` library with both features
enabled against the bounded original material source closure described above.

The Euler field remains frozen-composition. The zero-dimensional reactor does
not add species advection or chemistry/flow splitting to that field. Neither
path supplies diffusion flames, generic third-body/falloff mechanisms, spray,
knock, real-gas phases, moving volumes or a validated ignition model.

### Constant-pressure reactions on the common material parcel

`material::isobaric::MaterialIsobaricReactor` connects the existing prepared
constant-pressure chemistry owner to that same `ClosedGasParcel`. Its adiabatic,
held-bath and finite-body entry points reuse the original chemical stages,
temperature-dependent Cp, reaction formation enthalpy and entropy accounts.
The bath/body heat is integrated with the reactions; it is never an additional
chemical heat-release term or an independent temperature update.

The original material enthalpy replacement resolves the final products at the
original mass and pressure. Actual `Delta H = Q` and `Delta U + p Delta V = Q`
are checked along with the chemical endpoint's H, U, volume, work, temperature
and species fractions. Finite bodies retain their source phase plateau and
density; combined gas/body H and U-plus-both-volume-works close before either
endpoint publishes. Source/transport refusal or final cancellation preserves
the complete gas parcel and leaves the supplied body unchanged. Zero duration
preserves the parcel exactly after admission.

Kinetic preparation and dense feed storage are reused. The returned chemistry
view borrows that workspace; existing `to_owned` supplies an independent
snapshot when needed. Material endpoint resolution still allocates. The source
identity and neutral dilute-gas restrictions above remain. Absolute
cross-composition energy resolution must meet the requested Joule tolerances;
unresolvable reference offsets refuse. This is a homogeneous quasistatic
pressure constraint, not load inertia, spatial combustion or experimental
validation. The focused regressions in `material/tests/isobaric.rs` cover
independent Cp/Cv balances, an analytic temperature-dependent reaction/bath
trajectory, finite-body melting and both expansion works, exact replay,
source/domain refusal and early/late cancellation.

### Heated and cooled reactions in rigid material parcels

`material::isochoric::MaterialIsochoricReactor` connects the prepared rigid
chemical stages to the same source-bearing `ClosedGasParcel`. Adiabatic,
held-bath and finite-body calls reuse original temperature-dependent species
Cv and formation energy. Signed heat and composition evolve together; the
material owner reconstructs their one final endpoint with
`replace_model_with_heat`, without requiring a separately heated feed or an
adiabatic product intermediate. All final transport and acoustic properties
therefore use the accepted products, temperature and pressure.

Mass and volume retain their original bits. Actual material `Delta U = Q`,
initial/final chemical U, mass, density, temperature, pressure and generated
fractions must satisfy the declared cross-owner tolerances. A finite body's
original Cp/enthalpy chart, latent inventory and expansion work remain in its
existing owner: `Delta Ugas + Delta Hbody` and
`Delta Ugas + Delta Ubody + Wbody` are checked separately before publication.
Reaction formation energy is included once. Gas pressure work is zero; body
enthalpy is never relabelled as rigid-gas storage.

Source mismatch, unsupported product transport, exhausted work and final
cancellation leave the complete gas parcel and input body unchanged. Zero
duration skips material re-inversion after admission and preserves both
states exactly. The prepared chemistry and dense feed storage are reusable;
material endpoint resolution still allocates. Neutral dilute-gas source and
absolute energy-resolution restrictions above continue to apply.

The three regressions in `material/tests/isochoric.rs` cover independent
variable-Cv/variable-Cp reaction and heat trajectories with both implicit
methods, heating and cooling, body expansion work, updated material response,
unchanged inputs, repeatability, source/domain refusal and early/late
cancellation. Pinned read-only Rust parsing and scoped whitespace checks
passed; these new native regressions remain unexecuted because the shared
environment lacks sibling checkouts and available compilation memory.

### Caller-owned end ports and local cell sources

`Euler1d::step_with_ports(max, controls, ports, cancelled)` advances a
non-periodic interval whose end fluxes are owned by the caller. The SSP-RK2
predictor offers stage 0 at `t^n` and the corrector stage 1 at `t^n + dt`, each
as a `PortStage { stage, time, dt, left, right }` holding the end cells'
conserved and primitive states; the closure returns `[Option<FaceFlux>; 2]`
(`None` keeps the declared boundary). A returned flux is the whole end
exchange: the accepted report's boundary changes are `dt * area * flux`. A
refusing port refuses the step (subject to the declared halvings) and leaves
field and time unchanged; periodic intervals refuse ports.

The relations a port needs are public and source-owned:

- `FrozenGas::end_wave(interior, normal, pressure, cancelled)` returns the end
  state reached from `interior` along its outgoing wave curve (isentrope in
  expansion, Hugoniot in compression), its velocity derivative with respect to
  the end pressure and its signal speed. Its weak limit is `du = -n dp/(rho c)`
  and it stops gas at the wall solver's pressure.
- `FrozenGas::nozzle_exit(donor, back, cancelled)` is the steady isentropic
  exit of a stagnation donor into a back pressure, with choking; its
  incompressible limit is Bernoulli and its stagnation enthalpy is the donor's.
- `FrozenGas::stagnation_to_static_state(T0, pressure, velocity, cancelled)`
  recovers the static face state from `h(T0)-h(T)=velocity^2/2`, with signed
  axial velocity and caller-owned pressure. It reuses the reservoir's bounded
  inversion of original source enthalpy differences and temperature-dependent
  Cp. Density, sound speed and conservative-state precision are admitted by
  that same source. This closes energy after adiabatic jet mixing without
  imposing an isentrope or replacing the caller's pressure/choking law. An
  exhausted cold source path refuses. Analytic constant/variable-Cp tests
  cover both velocity signs, formation-reference invariance and cancellation.
- `stagnation_enthalpy(T)`, `face_flux(face, signal)` and
  `inflow_flux(face, T0, signal)` (energy flux = mass flux times the donor's
  stagnation enthalpy) close a face from one state.

The shared nozzle closure resolves the imposed pressure difference even near
equal reservoir pressures. It factors the temperature decrement by the original
source's isothermal entropy increment `R ln(p0/pb)` and evaluates the source's
Cp/T and enthalpy secants in that coordinate. The residual is relative to the
driving increment; an absolute Cp-scale tolerance cannot introduce a flow floor.
The derivative limit retains finite nozzle kinetic energy when the temperature
decrement is smaller than one representable absolute-temperature step. Returned
temperature and EOS state still have ordinary floating-point rounding.

Positive source Cv gives a conservative subsonic shortcut: when
`2 ln(p0/pb) < pb/p0` and the bounded cooling interval is source-supported, the
closure avoids the preliminary sonic solve. Stronger pressure drops retain the
existing source-owned sonic limit and choking decision. Neither route freezes
gamma or extends the temperature domain. The same closure serves public frozen
gas end ports and finite-volume network valves. Focused regressions compare
constant/variable-Cp analytic isentropes down to adjacent pressures, the exact
zero-drop limit, reversed donor flow, energy-reference invariance, cold-source
refusal and cancellation. A counter over the actual FrozenGas source compares
the complete weak solve with the previously mandatory preliminary sonic solve;
this is a source-query work comparison, not a throughput or scaling benchmark.

`Euler1d::apply_cell_sources(update, cancelled)` replaces every cell average
through a caller-owned local update between accepted steps (wall friction, wall
heat exchange). Every candidate is decoded through the original caloric law
before any is published; one refusal leaves the whole field unchanged. Time
does not advance.

Tests (`dynamics::evolve::port_tests`, 6): wall ports reproduce reflecting ends
bit for bit with ordered stages; an injecting port is the whole end exchange and
conserves mass to 1e-15 and energy to 1e-14 relative; refusal and periodic
refusal; the end wave's acoustic limit and wall pressure (1e-9); the nozzle's
Bernoulli limit (1e-4 at a 1e-4 pressure drop) and face-flux consistency with
the Riemann interface; transactional and exact cell sources. Measured on trj
lane t, log 20261030T190307Z (61 `dynamics` tests pass).

### Material cards in the shared molecular transport law

With both `material` and `experimental-euler` enabled,
`material::transport::MaterialTransport<L>` implements the existing
`multicomponent::dissipation::TransportLaw`. This connects the original
`ThermochemicalGasModel` transport cards to every consumer of that shared law,
including molecular species flow, reacting flow, finite-plenum networks and
the frozen Cartesian `fs-flux` consumer. No separate mixture formula or Euler
material registry is introduced.

Coefficient ownership is explicit: the original material evaluator supplies
viscosity and thermal conductivity using its source selection and Wilke/WMS
owners; the separately supplied `L` supplies only the complete canonical
binary-diffusivity matrix; bulk viscosity is an explicit nonnegative input.
The scalar coefficients returned by `L` are deliberately unused and never
serve as missing-card fallbacks. Pair transport retains its own source domains
and refusal behavior. Material cards do not imply binary diffusivity or bulk
viscosity. Admission compares the entire original NASA-9 species models, not
their labels, and requires the material bank to cover the complete pair axis.

The adapter retains the original material bank, cards, selection policies,
water/pack identities, complete external `QueryPoint`, and binary law. Each
actual face temperature and pressure replaces only those two coordinates in
the query, with canonical absolute-temperature/pressure meanings. Other axes
and their quantity conventions remain unchanged. Exactly normalized canonical
mole fractions pass through the original `with_composition` owner unchanged.
Native conservative faces may close a compensated sum while their canonical
ordered sum differs from one by roundoff. Only the native `64*epsilon*n`
unit-sum envelope is admitted; the original `with_mole_amounts` generated-simplex
owner reconciles those values. The actual resulting canonical fractions are
then supplied to BOTH material and pair-transport owners. Genuinely nonunit
inputs refuse. Exact-zero species are omitted from material queries; no trace
floor makes an absent species appear. The full pair matrix remains present because the
original diffusion owner supports evolution into those species. A newly active
species must have its real transport card support at the actual state.

`resolve` exposes the original resolved material state and query receipts for
the same face inputs. Source dimension/state admission precedes both material
and diffusion queries. The caller's binary output is copied only after complete
material resolution, successful original pair evaluation, shared matrix
validation and a final cancellation poll. Failed or cancelled evaluations leave
that output unchanged. Material refusals remain structured as
`DissipationError::Material(MaterialError)`, including through the existing
reacting-flow transaction. The adapter contains no accepted fluid state.

The original material resolver allocates and evaluates source cards at each
call; this path does not inherit the allocation-free claim of `KineticTransport`.
Five focused G0/G1/G3/G4 tests exercise complete query/source retention, zero
species, coefficient ownership, independent Fourier/viscous mode responses to
data-only card changes, no-query input admission, late output cancellation and
binary refusal, canonical three-active-species/one-zero roundoff reconciliation
through the same material and binary composition, and whole reacting-flow rollback when a new product's card
cannot support the final temperature. Fixtures are manufactured engineering
data, not measured gas validation. Existing source domains, neutral ideal-gas
assumptions, diffusion omissions and each consuming solver's numerical limits
remain in force.

`MaterialBinaryTransport` supplies the complete molecular law directly from the
original material gas and explicit `BinaryDiffusionMaterial` pair assignments.
The original caloric source is compared in full; each unordered species pair
must occur exactly once, including currently absent species. Pair cards retain
their original external contexts and pin the selected binary-diffusivity claims
at their declared reference T,p. Actual face T,p refresh those same claims;
missing support refuses rather than switching to another source or estimating
D from viscosity. Binary data have units m²/s and must be positive. Pair contexts
cannot freeze mixture-composition coordinates or aliases of evolving T,p.
Source axes and applicability bounds require supported static external meanings;
untyped restrictions and frequency/field/rate support refuse. Unused context
labels remain in the original query receipts, and declared external coordinates
such as humidity remain fixed while T,p evolve.

This law enters the existing axial, planar, spatial, finite-parcel and reacting
owners through `TransportLaw`; it does not introduce another transport solver.
The shared continuum owner retains its mass frame, pressure-driven separation,
original partial-enthalpy flux and full energy accounting. The original
`MaterialSoretDufourTransport` may wrap this law with independent cross data.
Scalar viscosity/conductivity still come from the original material gas, and
bulk viscosity remains explicit. Failed evaluations preserve the output and
last complete pair receipts. Source queries allocate and carry their existing
validity limits; this path is not an empirical validation or performance claim.
G0/G1/G3/G4 tests in `material::tests::transport::binary` check unequal-mass
composition/pressure diffusion and partial enthalpy against independent binary
answers, conservative spatial mode refinement, data-only diffusivity changes,
pair symmetry/coverage, source pinning and cancellation.

### Material bulk viscosity in compressing gas

`MaterialBulkViscosityTransport<L>` wraps any original `TransportLaw` with an
explicitly assigned mixture bulk-viscosity card. The source supplies the local
Newtonian coefficient zeta [Pa s], including a supported exact zero. Its
nonnegative value replaces only the wrapped law's bulk viscosity; shear
viscosity, conductivity, binary diffusion and Soret/Dufour heats retain their
original owners. Existing gas consumers use this coefficient in compressive
stress and its conservative energy flux. No extra viscous heating source is
added.

Every canonical species has an explicitly typed mole-fraction coordinate.
The reference absolute T,p selects and pins the original source claim; actual
T,p and the complete canonical composition refresh that same claim, retaining
external conditions and original receipts. Shared simplex admission permits
only native unit-sum roundoff. Source exits, negative coefficients, wrong
units, coordinate aliases and nonstatic support refuse. A frequency-dependent
apparent viscosity cannot replace relaxation memory. Callers composing a
separate relaxation channel must ensure the coefficient does not include the
same loss twice. Pure-species data do not imply a mixture law, and matching
dimensions do not make shear viscosity a bulk coefficient.

Identical-state queries reuse the last successful bulk resolution. Material
queries may allocate. Refusal/cancellation preserves caller outputs and the
last complete bulk receipt; wrapped workspaces/receipts are not accepted gas
state. G0/G1/G3/G4 tests in `material::tests::transport::bulk` exercise live
temperature/pressure/composition sources, compression and trace-free stresses
with conservative work, bulk-only kinetic-to-thermal conversion at constant
total energy, exact zero, coefficient/heat delegation, source admission and
atomic refusal/cancellation. Manufactured fixtures verify coupling and
numerical behavior, not empirical mixture calibration or a relaxation model.

## Reusable finite-rate reactor for cell-wise chemistry

`chemistry::IsochoricReactor` owns one original immutable kinetic mechanism and
one declared set of numerical controls. Construction admits the complete
caloric source interval and allocates the original source, stoichiometric,
Newton, adaptive-stage, concentration and endpoint workspace once. Subsequent
`advance_adiabatic` and `advance_adiabatic_energy` calls take dense `Amount`
values in canonical source order and reuse that storage without successful-path
allocation. This removes repeated identity lookup, source-path admission and
workspace construction when processing a finite-volume field.

The returned `IsochoricKineticsView` borrows the accepted endpoint and its
source owner. Explicit `to_owned(cancelled)` alone constructs the independent
result and original NASA-9 evaluation receipts. The existing allocating free
functions delegate to the same evolution algorithm. Each new solve resets
extents, extensive inputs and statistics; cancellation, an unsupported energy
target or a work-limit refusal cannot contaminate the next cell.

The previously declared implicit-Euler/step-doubling accuracy, positivity,
energy, element, mass and entropy boundaries remain in force. Focused native
verification passes all eight existing reactor regressions and four prepared
workspace tests covering analytic independent parcels, retained buffer
addresses, explicit-energy/source identity and recovery after refusals. This
is not a full-workspace or DSR verdict.

## Conservative transport of evolving ideal-gas species

Behind `experimental-euler`, `multicomponent::MixtureGas` owns one complete
original `IdealGasEquilibriumV1` species bank and an explicitly admitted
positive-Cv caloric interval. Every cell evolves all species partial mass
densities, axial momentum and total energy. Density is the compensated species
sum; temperature, pressure, Cv and frozen sound speed use the current
composition and original NASA-9 laws. Zero species stay admissible without
clipping, seeding or normalization.

The stored linear energy gauge is `rho E - sum(rho_i e_i(Tmin))`. Shared species
fluxes carry their corresponding original energy references, preserving total
physical energy while avoiding loss of sensible heat to large formation
constants. A two-species regression confirms bit-identical numerical transport
under independent formation-reference shifts of magnitude 1e18.

`MixtureEuler1d` uses the stationary `Mesh1d` owner, including nonuniform spacing,
face areas and integrated volumes. Its MUSCL/HLLC, bounded HLL fallback, SSP RK2,
actual stage CFL, pressure-area force, periodic/transmissive/reflecting boundaries
and transactional retries follow the existing engineering flow lane. Species,
momentum and energy use common face wave speeds. `MixtureGas::hllc_flux` exposes
the same source-bound Riemann closure with caller-supplied output storage.
Contiguous state, stage, slope and face workspaces are allocated once; successful
numerical steps allocate no workspace. Source or numerical refusal and
cancellation publish neither a partial field nor a new clock.

Six focused native regressions cover variable-Cp mixture thermodynamics and
sound speed; stationary differing-caloric contacts on a tapered nonuniform mesh;
smooth species-contact refinement; agreement with the frozen-source shock
solver at uniform composition; formation-reference invariance; and rarefaction,
structural zeros and late cancellation. Smooth contact L1 errors at 40/80/160
cells are 5.266259985e-4 / 1.198517512e-4 / 2.972969663e-5. The complete focused
pure-gas library passed 33 tests when this transport increment was admitted.

This remains an Estimated inviscid engineering lane. Fully conservative moving
material interfaces with differing caloric laws can exhibit pressure
oscillations; exact moving-interface pressure-equilibrium preservation is not
claimed. There is no fully discrete entropy-stability certificate, physical
species diffusion, viscosity, general moving geometry, arbitrary junction or multi-D flow
claim.

### Shared planar faces for evolving species

`MixtureGas` also admits borrowed `NormalCell` and `NormalPrimitive` traces.
The canonical output axis is every original species, normal momentum,
transverse momentum, then full species-referenced energy. Its axial and planar
consumers share one HLLC/HLL implementation; transverse momentum is not appended
to an independently solved one-dimensional flux. Both kinetic components enter
star-state and HLL intermediate-state admission through the original caloric
owner. In particular, shear mixing cannot hide heat outside the source interval.

The reconstructed-primitive entry reuses thermodynamics computed by the original
source constructor, avoiding an immediate second temperature inversion. All
successful face evaluations use bounded stack scratch and caller output storage.
Stationary slip walls retain the original composition-local shock/rarefaction
pressure; species, tangential momentum and energy wall fluxes are zero. Output
is scratch until success; neither input trace nor the species source is mutated.

Five focused native regressions cover exact axial parity, analytical material
contact transport, transverse-boost covariance on HLLC and HLL, formation-gauge
invariance, source refusal under shear heating, primitive/conservative parity,
wall pressure and cancellation. These are source-kernel claims, not a spatial
convergence, entropy-stability, chemistry or physical diffusion theorem.

## Coupled finite-rate multicomponent flow

`multicomponent::reactive::ReactingEuler1d` binds one conservative mixture field
to an identical original species bank in `IdealGasKineticsV1`. It composes an
isochoric chemistry half step, a complete conservative transport step, and a
second chemistry half step. If heat release changes acoustic speed enough for
transport to shorten the trial duration, the complete composition restarts
from the original field using the actual shorter duration. Failed chemistry,
source reconstruction, work budgets and cancellation publish no cell or clock.

Each chemistry solve uses the actual physical cell volume and total internal
energy after removing kinetic energy. Momentum and physical total energy are
held fixed. The stored energy changes only by the negative species-reference
change: `Delta Eg = -sum(e_i(Tmin) Delta rho_i)`. Original formation energy is
therefore accounted for once; no additional heat-release term is added.
The accepted reactor temperature is independently checked against the original
species sensible-energy functions at the reconstructed cell composition. Lost
absolute-energy precision, molecular-mass inconsistency that changes kinetic
energy, or an endpoint outside the source/energy budget is refused.

Two complete preallocated flow instances provide accepted and candidate state;
one prepared reactor and dense amount buffer are reused across cells. Source
identity/path admission and numerical workspace allocation occur at
construction. Chemistry tolerances are per physical cell; callers must choose
absolute amount and energy tolerances appropriate to the actual mesh volumes.
`ReactiveStepReport` retains transport boundary changes, area-wall impulse,
actual chemistry species/gauge-energy changes, chemistry entropy production,
solver work and mass/chemical-total-energy residuals. Chemistry entropy alone
is not a total transport/source/boundary entropy certificate.

Four coupled regressions cover a moving exothermic periodic field against the
homogeneous reactor at the actual accepted duration (including forced acoustic
restarts), independent sinusoidal advection/reaction decay and refinement,
late cancellation followed by reuse, and source-identity/chemical-energy
precision refusal. Reactive-contact L1 errors at 40 and 80 cells are
2.3273859838e-4 and 5.1358044353e-5. The complete current pure-gas native library
passes 38 tests. The five existing material-adapter regressions also pass using
actual fs-gas source with the selected original fs-material source bodies;
this is not a full-workspace or DSR verdict.

The runnable `reacting_oxygen` example uses the O/O2 high-temperature NASA7
coefficients and isolated `2 O + M <=> O2 + M` reaction from GRI-Mech 3.0 as
distributed in Cantera v3.2.0 `data/gri30.yaml`. NASA7 embeds exactly in NASA9
with zero inverse-temperature Cp terms. It retains the original 1000..3500 K
region, 1 atm standard-state pressure, explicit SI rate conversion and unit
O/O2 collider efficiencies. This isolated oxygen subsystem is not a complete
air or fuel-combustion mechanism and carries no experimental-validation claim.
A 40-cell periodic run reaches 1 ms in 113 accepted steps, with reported mass
residual 0 kg and chemical-total-energy residual -1.718625e-13 J; atomic oxygen
falls from 4.237578e-5 to 3.380236e-5 kg. Local temperatures reach
1831.691..1912.266 K from an initially uniform 1800 K.

Run with `cargo run -p fs-gas --features experimental-euler --example
reacting_oxygen -- 80`. The example emits coherent-SI CSV and a balance summary.
The symmetric composition does not upgrade the underlying adaptive
implicit-Euler chemistry to a certified second-order time method. Local
chemistry error estimates are not global coupled error bounds. Physical
species diffusion, heat conduction, viscosity, flames, sprays, real-gas phases,
general moving meshes and arbitrary reacting junction networks remain separate work.

### Finite-rate chemistry driving a resolved inertial piston

`ReactingEuler1d::step_with_piston` composes the original cell chemistry half
steps around `MixtureEuler1d::step_with_piston`. The first source uses the old
cell volumes; the second uses their actual moved volumes. Changes in species
and temperature therefore drive the original pressure waves and inertial wall.
The reactor's existing physical-energy target and species-reference correction
supply chemical conversion once; no additional heat-release source is added.
The constant-area interval has a fixed left wall and one inertial right wall.

Every complete attempt starts from the accepted field, geometry and piston
velocity. A shorter mechanical/acoustic duration restarts both chemistry half
steps at that duration. Second-half source refusal, failed whole-step balance
or cancellation cannot publish the successful middle pressure/motion stage.
The existing two field workspaces, reactor and amount buffer are reused. A
bounded, cancellable same-owner geometry reset copies only private trial
coordinates/volumes after motion; the ordinary stationary scratch-copy and
molecular-owner geometry guards remain intact.

`ReactivePistonStep` retains the original reactive and piston reports and adds
actual physical chemical-gas energy change plus whole-step momentum/energy
residuals. Endpoint species-reference changes participate in the physical
energy balance, together with gas energy, piston kinetic energy and exterior
pressure/additional-force work. The fixed left support carries its explicit
opposite impulse with zero work. The same caller-authored piston tolerances
gate whole-step balances as well as the pressure solve; callers must account
for their per-cell chemistry tolerances when choosing that global budget.

Four focused G1/G4 regressions cover unforced exothermic pressure-driven
expansion and species/work/support balances; acoustic restart with independent
stage replay at the accepted duration and moved volumes; a second chemistry
refusal after successful compression; and cancellation/recovery after prior
motion. They are authored but unexecuted under the session's worker resource
limit. This entry does not advance the separately bound fixed-geometry
molecular, finite-wall heat or network owners. Original chemistry method,
splitting-error, source precision and unequal-caloric-contact limits remain;
there is no general second-order, entropy-stability or combustion-validity claim.

### Molecular stress and heat transport through a reacting inertial piston

`ReactingEuler1d::step_with_piston_and_dissipation` evaluates the original
`D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2)` composition inside the shared
midpoint-speed piston root. Every proposal restarts all five stages from the
same accepted field at one duration and wall speed. The molecular stages use
the old and moved physical widths/volumes, respectively; only the original
extensive Euler A stage changes geometry and time. Chemical conversion,
Maxwell--Stefan species enthalpies, Fourier/Soret/Dufour transport and the
longitudinal Newtonian stress remain with their existing source owners.

Both pressure and molecular wall stresses drive piston inertia. Molecular
wall gradients use the candidate wall velocity, and the very same signed
momentum flux supplies laboratory wall work `w*(-tau)` in total energy.
`DissipationReport` now retains separate signed left/right conservative
transfers, including network branch reports. The returned
`PistonStep<ReactiveStepReport>` combines those impulses/work with the moving
Euler pressure exchanges. Actual endpoint chemical gas energy, piston kinetic
energy, fixed-left support impulse and exterior-pressure/additional-force work
must satisfy the caller's original physical budgets before publication. No
separate heat-release or viscous-heating source and no energy correction is added.

The constant-area interval has an impermeable adiabatic fixed left wall and
one impermeable adiabatic inertial right wall. Conduction and species diffusion
act inside the gas; there is no conductive wall-reservoir heat to count as
mechanical work. Finite wall storage, prescribed thermal end reservoirs,
variable-area or multidimensional motion remain outside this entry.

The molecular operator retains its accepted geometry separately from one
fallibly allocated candidate mesh. Trial widths are refreshed for every
stage, so a refused moved candidate cannot contaminate a later stationary
step. A shortened acoustic duration restarts the complete composition and
root. Gas, mesh, clock, piston and the molecular geometry binding publish
together after final acceptance and cancellation polling. Source and molecular
scratch are reused; no complete field is cloned per velocity proposal. The
original stationary molecular entry still requires its exact accepted mesh.

Four authored G1/G3/G4 regressions cover an independent one-cell RK2 viscous
wall force/work oracle; reacting expansion with nonzero viscous load, physical
energy/support balance and source-gauge changes; complete replay at the solved
speed after acoustic shortening; and final cancellation after prior motion,
scratch reuse and ordinary stationary continuation on the accepted geometry.
They are unexecuted in this resource-limited session. Original parabolic,
source, chemistry/splitting and unequal-caloric-contact limits remain; no
experimental combustion, general second-order, entropy or FSI-stability claim
is made.

### Finite end-wall heat driving reacting piston pressure

`ReactingEuler1d::step_with_piston_dissipation_and_heat` adds the original
caloric heat exchange around the same molecular/reacting midpoint-speed root:
`H(dt/2), D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2), H(dt/2)`.
`PistonThermalEnds` retains two distinct finite wall stores or prescribed
thermal reservoirs, their original caloric bindings and accepted gas geometry.
Every velocity proposal resets both stores from accepted values; its first
heat half uses old end-cell volumes and its last uses the moved volumes.
Heat therefore changes the actual caloric state used by pressure fluxes and
continuing motion. End order is left then right, including when both ends
contact the same gas cell. Only A advances geometry and time.

Pressure and longitudinal viscous stress remain the only gas-side mechanical
impulse/work terms. Conductive heat never enters the stress-work ledger. The
independent complete energy admission is
`Delta E_gas + Delta U_finite_walls + Delta K_piston + P_ext*A*Delta x
- F_ext*Delta x - Q_prescribed_reservoir = 0`.
Gas chemical energy includes actual endpoint species-reference changes once.
Finite-wall storage is reconstructed from accepted endpoints using the original
energy difference or nonlinear capacity integral; a cached capacity tangent is
not multiplied by the whole temperature change. Arbitrary wall energy datums
cancel in differences. Prescribed reservoir heat is external, whereas heat
exchanged with a finite store remains internal to this energy balance.

`ThermalReactiveStepReport` retains complete gas endpoints and heat transfers
alongside the original middle D/C/A/C/D report. That middle report's balances
exclude the outer heat halves. `PistonStep::total_energy_residual_j` includes
the finite stores and subtracts external reservoir heat on this entry;
`mechanical_energy_residual_j` continues to contain mechanical work alone.
The end owner's last report, contacts and accepted geometry publish atomically
with gas, clock, piston and molecular geometry after the final cancellation
poll. Acoustic shortening restarts all seven stages and the whole root.
Two small contact arrays and a prepared candidate mesh are reused; no thermal
source tables or complete fields are cloned per velocity proposal.

Source-backed storage must have an explicit internal-energy or admitted
negligible-expansion solid convention. Bare source enthalpy refuses because
its pressure/deformation work has no owner in this rigid-body balance. Radial
drag also refuses, avoiding a duplicate or stationary-wall momentum sink on
the moving piston. End conductance is caller-prescribed hA; finite thermal
masses represent distinct fixed amounts and do not add a second inertial mass.
Thermal expansion, axial wall conduction, latent material work not represented
by the admitted source, changing area and a fully implicit conjugate thermal
solution are outside this entry. Backward-Euler heat halves do not establish
global second-order accuracy or an entropy/combustion-validity theorem.

Four authored regressions cover an independent two-body heat oracle and hot/
cold finite-wall pressure response with full storage/work balance; nonlinear
original wall storage with simultaneous chemical and wall energy-gauge shifts
and a separately counted external reservoir; final cancellation after prior
motion, exact replay and original source-interval refusal; and admission of
the rigid thermal-work/drag convention. These tests are unexecuted under this
session's worker resource limit; no native build or solver result is claimed.

### Two resolved gas chambers sharing an inertial piston

`MixtureEuler1d::step_with_shared_piston` couples the original evolving
multicomponent fields on BOTH sides of one impermeable inertial wall. Each
chamber retains its own source, canonical species axis, caloric domain and
resolution. Their constant cross sections, shared interface coordinates and
physical clocks must match. The wall has negligible axial thickness; both
outer ends are fixed and reflecting. Positive axial motion expands the left
volume and compresses the right volume.

Each proposal in the original midpoint-speed root evaluates both original
moving-volume RK2 fields at one common duration from their accepted states.
The pressure on each side is its resolved moving-wall wave pressure, never
an imposed opposing reservoir or a guessed whole-chamber pressure law.
Piston impulse is the opposite of the summed interface gas impulses. The
very same face exchanges supply gas work; total gas energy plus piston
kinetic change equals the applied additional-force work. Fixed support
reaction is the opposite of BOTH outer gas impulses, with zero support work.
`SharedPistonGasReport` retains the two separate original conservation and
boundary reports, including their independent species balances and gauges.
The outer `PistonStep` sums fixed support reactions; exterior-pressure work
is zero because both gas sides are already represented.

The initial timestep obeys both relative acoustic limits and both possible
mesh-contraction limits over the authored speed bracket. A refused stage or
root halves and restarts the complete two-field proposal. Each original
field reuses its candidate/RK/geometry workspace; no source tables or full
fields are cloned per velocity proposal. After the last cancellation poll,
both extensive fields, meshes, equal clocks and piston velocity publish
together. Shared candidate endpoints must represent exactly one position.
`gas_trials` counts paired proposals, each requiring two original field solves.

Four authored G1/G3/G4 regressions cover equal-pressure rest with different
gases and independent opposing p*A*dt impulses; unequal-pressure waves in
both directions, source-energy gauge invariance and full work/support balance;
fine-side CFL shortening with replay of both prescribed moving fields; and
late cancellation, continuation plus incompatible geometry/clock refusals.
Native compilation, formatting and runtime execution remain unverified after
the delivery executor became unavailable. This is sealed adiabatic Euler
transport, without advancement of separately bound chemistry, molecular or
thermal owners. Original unequal-caloric contact and boundary-order limits
remain; there is no general FSI stability or multidimensional-rung claim.

## Characteristic compact-reservoir boundaries

Frozen Euler `Boundary::Reservoir` retains a stationary conservative state in
the original gas owner's energy gauge. Its pressure and temperature supply
stagnation conditions for inflow; its pressure supplies the outflow backpressure.
A one-sided classical-wave compatibility relation determines the actual port.
A stationary two-sided ghost would add an artificial acoustic impedance and is
not used for this compact-volume model.

Incoming gas accelerates isentropically using the source's temperature-dependent
enthalpy and entropy. A safeguarded velocity solve matches the outgoing duct
wave; an admitted sonic state limits choked inflow. The energy flux carries the
reservoir's original total enthalpy exactly once. Outgoing gas follows the
thermally perfect shock Hugoniot or rarefaction invariant and carries its
resolved total enthalpy. The fixed interface samples departing supersonic waves
and sonic fans rather than overprescribing an outlet.

Scaled pressure-wave roots retain weak-shock/rarefaction increments. All caloric
differences come from the original source, without constant-gamma replacement.
Actual states and classical wave ordering must be admitted; cancellation and
bounded root work remain explicit. Unused stopped shocks or downstream fan
tails may lie outside the caloric interval without disqualifying an admitted
accelerating inlet or sonic outlet. A genuinely unsupported port state refuses.

`set_boundaries` updates prescribed end conditions atomically between steps.
Reservoir conditions are held fixed during each SSPRK2 step. The accepted
`left_boundary_change` and `right_boundary_change` are each signed into the
interval and retain the two-stage mass, momentum and energy exchange, alongside
the existing aggregate report. This supports driven inlets/outlets without
inferring exchange from a sampled pressure or flow value.

Seven focused native boundary regressions compare gamma-law shock/rarefaction
jumps and pressure derivatives (including 1e-10 pressure increments), analytic
sonic inlet/outlet limits, independent variable-Cp Hugoniot and quadrature,
reference-gauge invariance, source-window edge cases and driven-field balances.
All 28 frozen-dynamics tests passed with the integrated source. A separate
regression verifies that multicomponent flow refuses this frozen-state reservoir:
it requires a species-bearing state and energy-reference closure and must not
silently reinterpret this boundary as transmissive.

This estimated inviscid model assumes an ideal compact reservoir/entrance.
It supplies neither an empirical discharge coefficient nor finite neck inertia,
wall loss, multidimensional mixing, moving-volume work or an entropy theorem.

## Synchronized finite-plenum gas networks

`dynamics::network::GasNetwork1d` owns existing Euler branches, finite rigid
`Plenum` volumes, and explicit `PlenumConnection` endpoints. Every branch
retains its original mesh and field. All branches must share the exact frozen
source, composition, caloric interval, energy reference and accepted clock.
Duplicate, out-of-range and periodic endpoint connections refuse. Branch and
plenum counts are each capped at 4096; combined cell/plenum storage shares the
existing `MAX_EULER_CELLS` cap.

Each connected port evaluates the same physical characteristic reservoir
closure once at each stage. Its area-weighted mass and total-enthalpy exchange
feeds the branch and plenum with opposite signs. A compact plenum stores mass
and internal energy, without resolved bulk momentum; incoming kinetic energy
remains in that energy account. Its fixed support receives the reported port
impulse in the branch's axial coordinate and does no mechanical work.
Different branch axes are not added into a fictitious scalar network momentum.

One global SSPRK2 duration advances every owner. Actual branch face speeds and
the summed port area/speed divided by plenum volume limit both stages. A small
plenum can therefore limit the whole network timestep. Every retry restarts all
owners; state, clock and accepted exchange/report slices commit together only
after the final admission and cancellation check. Numerical scratch and report
storage are allocated at construction and reused.

Global reports integrate all branch and plenum mass/energy and count only
external-end exchange. Accepted per-port exchanges are the exact negation of
the corresponding branch-end report. Each branch also retains its full
momentum/area-force balance. Immutable access prevents independent branch
stepping while the network owns the common state.

Ten focused native network regressions cover original isolated-Euler identity,
tapered rest, compact-volume CFL, topology/source/clock admission, final-stage
cancellation and failure rollback, finite-pulse splitting, direct extensive
storage and port exchange reconstruction, orientation symmetry and bitwise
variable-Cp formation-gauge invariance. A 64/128-cell compact-plenum acoustic
check compares against `p'(t)/p'(0)=exp(-c*A*t/V)` before returning reflections;
normalized errors are 2.97668e-6 and 3.61846e-6, both within the fixture's
weak-wave tolerance. This is an independent impedance check, not an assertion
of an asymptotic refinement order.

The `fs-duct` `nonlinear_manifold` consumer uses original cylinder/cone
geometry and a 30 mL plenum. Runs with 16 and 64 cells per branch reach 4 ms
in 131 and 527 accepted steps. Both report zero mass residual; physical energy
residuals are -1.421085e-14 J and zero respectively. All temperatures remain in
the declared N2 200--1000 K source window, and both daughter passages respond.
The full staged pure-gas native run passed 58 tests and the existing eight
focused duct tests passed. This is not a workspace Cargo/DSR result.

A finite plenum is an explicitly compact, perfectly mixed storage model, not
a zero-volume lossless junction or a resolved three-dimensional shock/mixing
region. Species-changing/reacting junctions, valves/discharge coefficients,
wall heat, viscous loss, compliant volumes and entropy certificates remain
separate models. This contributes to the open broad gas-network bead without
claiming its complete reactive/MoC scope.

## Physical reflecting walls for evolving mixtures

Reflecting multicomponent boundaries solve the one-sided stopped shock or
isentropic rarefaction with the incident trace's current species composition.
They share the thermally perfect wall kernel with frozen Euler instead of
retaining the momentum component of a mirrored HLLC/HLL ghost flux. Every
species mass flux and the stationary adiabatic wall energy flux are exactly
zero; the physical wave speed enters the existing stage CFL calculation.

The borrowed local caloric law uses the admitted original species models,
temperature-dependent Cp/Cv, molar masses, energy/entropy differences and
source seams. Composition is fixed across each local inviscid wall wave but
may change between faces and accepted steps, including through the existing
reactive transport owner. Successful wall evaluation allocates no heap storage
and constructs no replacement gas/source owner. Source failures remain typed;
unsupported hot/cold wave paths and numerical budgets refuse without
extrapolation or clipping.

Five focused regressions compare arbitrary constant-Cp compositions against
analytic shock/rarefaction pressures, variable-Cp compression against independent
mass/momentum/energy jumps and entropy production, and expansion against a
separate Simpson integral. They also check orientation and independent species
formation-reference invariance, actual finite-volume wall wiring, and evolving
species/energy conservation with wall impulse accounting. In the retreating
manufactured fixture, the previous mirrored flux gives -0.722497216 Pa traction;
the physical wall pressure is +0.0110814482 Pa.

The complete integrated pure-gas native run passes 67 tests, including existing
reacting flow, characteristic reservoirs and finite-plenum networks. Eight
existing nonlinear duct tests pass against the refactored kernel, and the
64-cell-per-branch manifold still reaches 4 ms in 527 steps with zero reported
mass and physical-energy residuals. These are focused production-code checks,
not a workspace Cargo/DSR result or experimental validation. Reactive/catalytic
surface kinetics and wall heat transfer remain separate models.

### Resolved moving-wall pressure and mechanical work

Frozen and evolving-mixture planar/spatial face owners additionally expose
`moving_wall_normal_flux` and `moving_wall_spatial_flux`, each with a
source-admitted primitive fast path. These are moving-face ALE fluxes. The
original one-sided shock/rarefaction owner uses normal velocity relative to the
wall and retains the same density, temperature, pressure, actual species
composition, temperature-dependent calorics and source-domain refusals.

In the supplied normal axis the flux contains exactly zero relative mass and
species transfer, normal momentum `p_star`, zero tangential momentum, and
laboratory-frame energy `p_star * wall_velocity`. Multiplication by the existing
-1/+1 outward-normal sign gives outward transfer. Thus a piston doing work on
the gas adds internal and kinetic energy through its pressure traction, with
no duplicated heat or species-enthalpy transfer. Tangential translation has no
effect on the inviscid slip pressure or work. Independent species formation
references cancel because no species crosses the wall.

The reported signal-speed bound is relative to the moving face, not the
laboratory frame. Zero wall velocity retains stationary arithmetic and exactly
zero work. Both adapters use the existing shared caloric wave kernel; successful
face work creates no new gas/source owner and needs only bounded stack scratch.

Seven added regressions cover independent constant-Cp piston shock/rarefaction
pressures, variable-Cp swept-slab mass/momentum/energy balances, evolving-mixture
species and formation-reference invariance, normal/tangential Galilean changes,
relative wave-speed semantics, planar/spatial and conservative/primitive
agreement, stationary identity, source limits and cancellation without partial
work publication. They have not been executed in the current resource-limited
session. This remains an instantaneous inviscid face closure: moving geometry,
cell-volume evolution/geometric conservation, temporal quadrature, piston
inertia and viscous wall work are the corresponding spatial/mechanical owners'
responsibility. The separate moving-interval entry below supplies that spatial
evolution for a bounded axial geometry; ordinary stationary updates are unchanged.

### Spatial piston dynamics and conservative moving volumes

`Euler1d::step_with_moving_walls` advances the original resolved Euler field in
a constant-area interval with two impermeable prescribed end walls. Their two
axial speeds are fixed over each accepted step; all interior faces move by
affine interpolation, retaining the original nonuniform cell fractions.
The existing MUSCL reconstruction, source admission and HLLC/HLL wave owner are
shared with stationary dynamics. An interior Riemann problem is solved in its
face's translating frame, then its flux is returned to laboratory conserved
variables as `[f_mass, f_momentum+w*f_mass,
f_energy+w*f_momentum+(w*w/2)*f_mass]`. Thermodynamic traces are reconstructed
through the original source in that frame, avoiding subtraction of large
laboratory kinetic energies. Wall pressure and `p_star*w` work use the existing
one-sided caloric shock/rarefaction closure.

Both RK2 stages advance extensive `Q=V*U`, with the predictor resolved at the
new volume and the corrected endpoint on that same new geometry. Cell volume
changes are exactly the discrete swept volume `dt*A*(w_right-w_left)` used by
the moving-face flux. Actual face-relative CFL is checked at both old and new
stage volumes, and contracting geometry also bounds the initial trial. The
existing bounded retry/source-refusal policy never clips conserved fields.
`StepReport` retains the same stage-integrated left/right pressure impulses
and wall work that changed the gas; its conservation residual subtracts those
boundary exchanges once. Geometry scratch is fallibly allocated at the first
nonzero-motion entry and reused on subsequent steps. Cancellation/refusal
preserves the accepted mesh, cells and clock together; zero speeds retain the
ordinary stationary path and its exact arithmetic.

Six focused regressions cover independent Reynolds-transport and Galilean flux
relations, zero-speed identity, a co-moving nonuniform pressure-contact field,
resolved variable-Cp piston compression and formation-gauge invariance, a
16/32/64-cell homologous-compression comparison against the independent smooth
isentropic Euler solution, and early/final cancellation or geometry/source
refusal. They are authored but unexecuted in the resource-limited session.
The original unreconstructed wall-cell traces can limit global boundary order;
second-order global convergence is not asserted. This lane does not evolve
piston inertia, heat transfer, composition, radial/variable-area deformation,
open ends, bound gravitational atmospheres or multidimensional moving meshes.
The face-frame conservation convention follows Springel (2010), §3.3,
<https://doi.org/10.1111/j.1365-2966.2009.15715.x>; it does not claim that paper's
Voronoi mesh algorithm or its full accuracy properties.

`MixtureEuler1d::step_with_moving_walls` supplies the same bounded axial motion
for the original evolving-composition field. It shares the geometry workspace,
affine face velocities and swept-volume preparation with the frozen owner,
while retaining its own compatible mass-fraction reconstruction, canonical
species storage, source calorics and shared HLLC/HLL waves. Every partial
density advances extensively on the same RK2 volumes. Face-frame species
fluxes are unchanged by the laboratory transformation; their compensated sum
is the mass flux used in both momentum and kinetic-energy transformations.
There is no independently rounded mass-flux field, composition normalization
or absent-species seeding. Each species keeps its own formation-energy datum.
All wall species fluxes are zero, and the returned original mixture report
accounts for the pressure impulse and mechanical work at both ends.

The shared geometry's cancellation is translated into the mixture's native
`Cancelled` refusal before retry handling. Accepted mesh, all species, energy
and clock remain atomic. Five focused regressions cover independent species
Reynolds-transport fluxes, stationary identity, co-moving unequal-caloric
contacts, nonuniform-composition piston work and independent species gauges,
single-species agreement with the original frozen owner, structural zeros,
and early/late geometry or source refusal. These are authored but unexecuted
under the same worker resource limit. This standalone entry does not advance
the separately bound stationary reaction, diffusion, heat or network owners.
Generally advecting unequal-caloric material contacts retain the original
conservative scheme's pressure-oscillation limitation; no exact general
pressure-equilibrium or entropy-stability claim follows from this extension.

### Simultaneous inertia of both resolved gas end walls

`Euler1d::step_with_inertial_walls` advances a closed constant-area field and
two finite `Piston` masses together. Both wall velocities use the increasing
axial coordinate. Each `InertialWallLoad` holds exterior pressure and an
additional signed force: exterior pressure pushes the left wall rightward and
the right wall leftward. Gas momentum remains in the original resolved field.

A bounded two-variable Newton solve selects both midpoint wall speeds.
Every trial runs the original moving-volume RK2 operator from the same
accepted state and duration. Each wall receives the negative of its actual
gas-side pressure impulse; the identical boundary work enters its kinetic
balance. Both wall momentum/work residuals and combined gas-plus-walls
momentum/energy are admitted independently, using represented mesh
displacements for exterior/load work. No fixed support, reduced wall mass,
uniform cavity pressure or posterior energy correction replaces an owner.

Both velocity intervals bound relative-wave and contraction CFL. Row-scaled
finite differences and bounded backtracking reuse the native candidate and
geometry buffers. An exhausted solve retries the complete interval under
the original gas halving budget. Only the last evaluated native candidate
can publish; gas, mesh, clock and both wall velocities commit after the
final cancellation check. No field or source table is cloned per trial.

Six regressions cover balanced outside pressure, symmetric vacuum recoil,
unequal walls with nonuniform pressure and signed loads, a common velocity
boost through the original first-order ALE route, the massive-left-wall
limit of the original one-piston owner, and cancellation/refusal rollback.
They are authored but unexecuted under the session's worker resource limit.
This frozen-composition adiabatic lane does not yet compose the two moving
walls with mixture chemistry, heat transfer, gravity, contact stops or
transverse dynamics; it establishes no general FSI stability theorem.

### Evolving mixtures with two finite end walls

`MixtureEuler1d::step_with_inertial_walls` now uses the same simultaneous
wall solve and physical acceptance with its original species-resolved ALE
operator. Local partial densities and original caloric sources determine
pressure at every stage. The complete native gas report retains canonical
species masses, both boundary exchanges and independently referenced total
energy. Each impermeable wall has zero flux of every species; absent species
remain absent. No nominal frozen composition replaces the evolving field.

`InertialWallsStep<R>` retains the frozen report as its public default. A
small internal exchange record passes only actual gas momentum/energy changes
and each boundary impulse/work to the shared mechanical equations; no species
or source state enters a duplicate wall model. Trials reuse the native mixture
candidate and geometry workspaces. Gas, all species, mesh, time and both walls
publish only after the common physical checks and final cancellation gate.
Original mixture failures retain their typed causes; cancellation and
allocation failures propagate without numerical retries.

Three focused regression groups cover unequal walls with heterogeneous
composition and independent formation-gauge shifts; agreement with the
original frozen owner for a fixed two-species composition; and early/late
cancellation, failed roots and unsupported-boundary rollback. These tests
are authored and formatting checked, but have not been compiled or executed
in this constrained environment. Chemistry, molecular diffusion, viscosity
and wall heat remain separate owners. The underlying conservative mixture
scheme retains its pressure-oscillation limitation at generally advecting
unequal-caloric contacts.

### Finite-rate chemistry with two inertial gas end walls

`ReactingEuler1d::step_with_inertial_walls` composes the original source
half steps around the simultaneous two-wall mixture pressure solve. Both
wall momenta remain finite and gas momentum stays spatially resolved. The
first chemistry half uses accepted cell volumes and the second uses actual
moved volumes. A shortened native interval restarts both source halves and
both walls from accepted state at that same duration.

`ReactiveInertialWallsStep` retains complete reactive accounting and the
original middle-stage wall report. Whole-system momentum subtracts both
signed exterior-pressure/additional-force impulses, with no fixed-support
term. Whole-system energy uses actual chemical gas-energy change, both wall
kinetic changes and represented exterior/additional-force work. Each changed
species reference participates once through the original coupled report;
no separate chemical heat or posterior energy correction is added.

All source, field and geometry workspaces are reused. Species, mesh, clock
and both wall velocities publish together after complete physical acceptance
and the final cancellation poll. Four authored G1/G3/G4 groups specify
unequal-wall reactive recoil with an independent energy-gauge comparison,
accepted-duration replay using moved source volumes, postcompression source
refusal and cancellation/retry after prior motion. They parsed with the
pinned formatter but have not been compiled or run in the constrained
editing environment. This chemistry-only entry does not advance molecular
transport or thermal walls. Existing chemistry/splitting and unequal-caloric-contact limits
remain; no general accuracy, entropy or combustion-validity claim is added.

### Molecular stress and finite heat at both inertial gas ends

`ReactingEuler1d::step_with_inertial_walls_and_dissipation` composes the
original D/C/A/C/D stages inside the simultaneous two-speed root.
`step_with_inertial_walls_dissipation_and_heat` adds the original finite
end-store or prescribed-reservoir H stages around that same composition.
Every proposal starts from accepted owners; second halves use actual moved
volumes. The identical candidate velocities supply pressure work, viscous
gradients and molecular stress work at BOTH finite walls.

Each wall receives its full pressure-plus-viscous impulse. Native acceptance
checks both impulse/work equations and combined gas-plus-walls momentum.
Chemical gas energy uses actual species-reference changes once. Total
energy adds both kinetic changes, original finite-store energy differences
and represented exterior/load work, then subtracts prescribed reservoir
heat. Conductive heat never enters mechanical work, and thermal storage
adds no duplicate translating mass. Original material caloric conventions
and source domains remain authoritative.

The existing two-speed solver retains its original frozen/mixture wrapper
and arithmetic; a private generic-error seam preserves source refusals.
Shortened acoustic intervals terminate every root proposal and restart
all stages at one common duration. No field/source clone is made per
proposal. Gas, geometry, clock, both velocities, molecular binding and
thermal stores publish together after physical acceptance and a final poll.

Four authored regressions cover unequal-wall viscous/chemical recoil and
source-gauge invariance, independent sequential finite-heat and complete
energy oracles, replay after acoustic shortening, and cancellation/retry
after prior motion with both finite stores and geometry bindings retained.
They parse with the pinned formatter; native compilation/execution was
unavailable in this resource-constrained session. Existing source, explicit
parabolic, chemistry/splitting and unequal-caloric-contact limitations apply;
this composition supplies no general accuracy, entropy or FSI theorem.

### Pressure-driven inertia of a resolved gas end wall

`Euler1d::step_with_piston` couples the existing spatial Euler field to
`dynamics::piston::Piston`, a finite inertial mass at the right endpoint. The
left wall remains fixed, both ends are impermeable, and area is constant.
The mesh alone owns wall position and area. Authored exterior pressure and
additional force act on the piston; neither a uniform chamber law nor a
second pressure model replaces the source-owned wave field. The existing
zero-dimensional `fs-couple::cavity_phs::piston` owner remains independent.

Each bounded scalar trial chooses midpoint velocity `w`, then runs the
original moving-volume RK2 operator from the accepted field at the SAME
duration. Its accepted right-wall pressure impulse `I` enters
`M*(v1-v0)=I+dt*(F-A*p_ext)`, with `v1=2*w-v0`. The identical two-stage wall
pressure supplies gas work. Consequently `Delta K=w*M*Delta v` pairs the
mechanical update with the discrete gas exchange; no subsequent energy
correction is applied. The report independently measures actual gas energy,
wall kinetic change, exterior-pressure and force work, and both mechanical
and combined momentum/energy defects. The fixed-left-wall reaction is
explicit in the combined momentum balance and does zero work. External work
uses represented mesh displacement, so a rounding discrepancy is admitted
only within the same caller-specified physical energy tolerance.

`PistonControl` requires a source-supported midpoint-speed bracket containing
the accepted velocity, positive momentum/energy tolerances, and a bounded
secant/bisection budget. Bracket speeds bound the initial relative-wave CFL;
both gas stages retain original source, moving-geometry and CFL admission.
Any retry halves and resets the complete coupled attempt. Only the last
successfully evaluated candidate can publish gas, geometry, clock and piston
velocity, after a final cancellation poll. Original gas workspace is reused;
geometry scratch is allocated fallibly on first use. There are no per-trial
field clones, and cancellation is never converted into a numerical retry.

Five focused G0/G1/G3/G4 regressions cover balanced pressure, variable-Cp
compression with source-energy-gauge invariance and global work/momentum
closure, exact-duration replay of the final solved motion through the
prescribed-wall owner, three mesh refinements against the independent
small-amplitude acoustic piston-impedance solution before reflections, and
early/late cancellation or source/bracket/geometry refusal. These tests are
authored but unexecuted under the session's worker resource limit. The
acoustic oracle is a linear small-amplitude limit; neither it nor residual
closure establishes experimental validity or a general FSI stability theorem.
Existing spatial-order and unequal-caloric-interface limits remain unchanged.
This entry does not supply piston collisions, friction, springs, crank
linkages, multiple moving inertial walls, heat or reactive/multicomponent
coupling. Context: Lefrançois and Boufflet, SIAM Review 52(4), 747–767 (2010),
<https://doi.org/10.1137/090758313>; the discrete work identity above is derived
directly for this implementation, not inherited from that paper.

`MixtureEuler1d::step_with_piston` now supplies the same inertial end wall for
the original evolving-composition field. It uses the identical `Piston`,
`PistonLoad`, `PistonControl`, scalar root and mechanical balance code. A small
source-neutral exchange carries actual gas momentum/energy changes, both
wall impulses, right-wall work and elapsed duration into that shared code.
`PistonStep<R>` retains the native gas report: the public default remains the
original frozen `StepReport`, while the mixture method returns its original
species-bearing report. No species, total-mass or energy-reference axis is
collapsed into a frozen surrogate.

Every mixture trial calls the existing extensive species ALE operator at the
common nonlinear duration. Actual partial densities, independent formation
references and source calorics determine pressure throughout the field.
Closed walls have zero flux of every species, and structural zeros remain
zero. The complete accepted field, geometry, clock and piston publish only
after shared momentum/energy acceptance and final cancellation polling.
Original mixture errors remain distinguishable as `PistonError::MixtureGas`;
native source/geometry cancellation is converted to `PistonError::Cancelled`
before retry handling. The prescribed-motion mixture path shares only its
original report extraction and trial visibility; its arithmetic is unchanged.

Four additional G0/G1/G4 regressions cover heterogeneous-composition inertial
compression with every species balance, independent source-gauge shifts and
global work/support exchange; single-species agreement with frozen inertia;
replay of the final solved midpoint speed through the original prescribed
mixture motion; and cancellation plus distinctly typed source, bracket and
geometry rollback. They are authored but unexecuted under the same worker
resource limit. This entry does not advance separately bound reaction,
diffusion, heat or network owners. The existing pressure-oscillation limit at
general moving unequal-caloric contacts remains explicit; shared mechanical
balance acceptance is not an entropy or pressure-equilibrium theorem.


## Conservative molecular, thermal and viscous field evolution

`multicomponent::dissipation::DissipationOperator<L>` binds the original gas,
stationary constant-area mesh, boundary conditions and explicit transport law.
`advance_for` completes an exact requested duration through bounded SSP RK2
parabolic substeps. It publishes the complete field atomically and leaves its
clock unchanged for the composing integrator. Nonuniform spacing is admitted
with the actual distances and physical cell volumes; variable-area viscous
geometry requires additional stress terms and is refused.

Every shared face carries barycentric Maxwell-Stefan mass fluxes Ji, momentum
flux `-tau`, and, without thermal cross coupling, gauged-energy flux
`sum_i (h_i/M_i - e_ref_i) Ji - lambda grad(T) - u*tau`.
The longitudinal stress is `tau=(4mu/3+mu_bulk) grad(u)`; bulk viscosity is an
explicit authored value, with zero selecting Stokes' hypothesis. The source
enthalpy difference is evaluated as `[Delta u_i(Tmin,T)+RT]/M_i`, avoiding
large-formation-reference subtraction. No independent viscous-heating source
is added to total energy. Viscous loss of kinetic energy therefore appears as
internal heating through the conservative update.

Periodic faces are shared bit-for-bit. Transmissive ends use zero diffusive
gradients. Stationary reflecting ends are impermeable and adiabatic, retain
viscous traction, and perform exactly zero wall energy work. Frozen reservoir
states have no species axis and are refused by the mixture boundary owner.
Parabolic estimates are checked at both stages; positivity, source support,
work limits and cancellation use bounded refusal/retry, with no clipping of
conserved states. Complete species, momentum, physical-energy and boundary
balances are returned.

`TransportLaw` accepts actual face T,p and composition. `ConstantTransport`
retains explicit coefficients and validity for authored approximations and
analytic cases. `KineticTransport` adapts the source-bound Lennard-Jones,
Chapman-Enskog/Wilke/Eucken model with reusable property storage. Matrix,
face, stage and candidate buffers are prepared once. Derived fraction/gradient
simplex identities use a checked roundoff-sized remainder on the largest
component; no conserved species value or physical flux is projected.

Seven focused native tests pass: Fourier heat-mode decay/refinement; analytic
binary diffusion; viscous velocity decay and kinetic-to-internal conversion;
unequal-mass partial-enthalpy transport and bit-identical evolution after a
1e18 formation-reference shift; nonuniform geometry; wall traction with zero
energy work; and cancellation, work limits, zero duration and reuse. Fourier
L1 errors at 40/80 cells are 1.91446548e-5 / 4.77807730e-6. This remains an
Estimated continuum operator: the parabolic estimate is not a spectral or
fully discrete entropy/invariant-domain proof. Explicit source-owned Soret/Dufour
coupling is described below; nonideal phases, Knudsen transport and variable-area
viscous stresses remain separate closures.

### Prescribed composition and temperature baths on the axial field

`DissipationOperator::advance_with_reservoirs_for` and
`advance_stabilized_with_reservoirs_for` accept optional stationary
species-bearing left/right reservoir cells for an entire molecular interval.
They extend the original SSP RK2 and RKL2 owners through the same spatial face
evaluation. `Some` overrides that end's molecular condition; `None` preserves
it. Two absent reservoirs execute the original entry point exactly, including
its cancellation sequence. The raw borrowed `Cell` representation follows
`MixtureEuler1d::step_with_reservoirs`: callers must match the owner's canonical
species order and per-species energy gauges because it carries no source
identity token. The original caloric source admits the complete bath state and
heat capacity before any trial. Periodic geometry, nonstationary baths and
invalid source states refuse without field or clock publication.

Boundary gradients use each actual end cell's center-to-boundary half-width,
including nonuniform spacing. The same constitutive kernel as planar/spatial
reservoir faces evaluates Maxwell--Stefan, Fourier, Soret/Dufour and viscous
transport from the fixed bath and the current stage's interior state. Species
energy uses the bath's original referenced enthalpies, while traction work
uses its zero velocity. No cell-center midpoint work or separate enthalpy
source is added. Barycentric diffusion preserves zero net mass transfer while
individual species, momentum and heat may cross an end. The existing signed
left/right and complete conservative balances integrate those fluxes with the
same SSP weights or RKL2 recurrence as the field.

Bath traces and caloric capacity use bounded stack storage; no per-face or
per-stage numerical workspace is added. Reservoirs remain fixed and undepleted
and are not retained after the call. Complete-duration atomic publication,
bounded work/cancellation and unchanged composing flow clock are preserved.
These methods advance molecular transport only: a composing integrator must
apply its matching advective pressure reservoirs separately. Reactive/network
reservoir composition and moving or finite-capacity molecular baths require
their respective coupling owners.

Five focused regressions exercise both time integrators: independent
Dirichlet binary eigenmode decay toward the linear steady profile and signed
species balances; steady Fourier transport across a nonuniform mesh with
physical half-cell gradients; unequal-species bath enthalpy and zero stationary
traction work under a 1e18 formation-reference shift; late cancellation,
work/source/geometry refusal and reuse; and exact absent-bath parity. They
establish the specified discrete field behavior, not an experimental transport
validation or a global entropy/invariant-domain theorem.

### Stabilized molecular diffusion intervals

`Rkl2Scheme` is the shared bounded Legendre coefficient owner for axial,
network, planar and volumetric consumers. It admits 2–256 stages and exposes
the stability extension, first weight and complete signed stage recurrence.
All consumers retain their original spatial divergence and physical gates;
the coefficient object supplies no substitute constitutive or boundary law.

`DissipationOperator::advance_stabilized_for` optionally advances the same
stationary constant-area field with second-order Runge–Kutta–Legendre (RKL2)
super-time-stepping. Every intermediate and endpoint evaluation uses the
original Maxwell–Stefan, partial-enthalpy, Fourier, stress and Soret/Dufour
face owner at its actual state. The linear negative-real stability extension
is `(s*s+s-2)/4` relative to the existing parabolic estimate. Original physical,
caloric and source admission still applies; failed stages halve the trial,
with no clipping of species, temperature or flux.

The three existing state arrays are reused, with one lazily reserved initial
face array and no per-stage allocation. Signed left/right boundary transfers
follow the same RKL2 recurrence as the conserved field, including its negative
weights. Complete-duration publication and the unchanged composing flow clock
retain the original transaction. Cancellation, exhausted attempts or exhausted
RHS calls leave the caller field unchanged.

`StabilizedDissipationControl` bounds macrostep duration, stage count (2–256),
retry halvings, accepted/rejected attempts and actual RHS calls. The report
retains the original physical balances and adds work, accepted stage counts,
and the parabolic number scaled by the stability extension. Rejected trials
and final source qualification consume the RHS budget. The caller chooses an
accuracy-appropriate maximum macrostep; this path has no temporal estimator.
The extended linear interval is not a nonlinear stability, SSP positivity,
entropy-decay or unconditional-stability proof. Moving geometry, chemistry
and hyperbolic composition remain their separate owners.

Eight focused native tests cover an independent Legendre heat-spectrum oracle,
physical-time refinement, binary diffusion, unequal-species enthalpy and
Soret/Dufour reference invariance, wall transfers, nonuniform cell volumes,
nonlinear rate growth with real rejected trials, and late refusal/reuse.
They pass with the original molecular/planar/spatial regression subset
(33 tests), including the current Maxwell–Stefan composition implementation.
For the tested 128-cell heat interval, RKL2 uses 19 original RHS calls versus
164 for SSP RK2, including its endpoint check. This is an operation count,
not a wall-clock speedup or measured-material validation claim.

### Explicit reciprocal thermal diffusion

`TransportLaw::heats_of_transport` optionally supplies canonical finite kinetic
heats of transport `H_i [J/kg]` at the actual common face T,p,composition.
`SoretDufourTransport<L>` supplies a constant authored approximation on an
explicit T/p domain and the full composition simplex, with source identity,
version and typed species coefficients matched to the underlying caloric bank.
The wrapped molecular law retains its own source/domain checks. No coefficient
is inferred from species names, molecular masses, heat capacities or collision
diameters. State-dependent or composition-restricted laws implement that same
source query directly. Existing molecular laws omit the extra closure by default.

With `material`, `MaterialSoretDufourTransport<L>` supplies this closure from
original material-card scalar, curve or two-coordinate surface claims. Each
canonical species is explicitly assigned a signed J/kg kinetic coefficient
key, card and selection policy, and an explicit typed mole-fraction query axis.
The caller declares the wrapped law's conductivity at zero diffusive species
flux; partial enthalpies cannot substitute for these kinetic coefficients.
Construction resolves a reference state, then pins the exact selected claims.
Every face refreshes the original sources at actual T,p and composition,
preserving fixed external conditions and full usage receipts. Original static
source admission rejects unowned frequency, field or rate response and
unsupported scalar constancy. Missing/ambiguous claims, type/domain/source
refusal and cancellation publish no partial coefficients or new retained cross
receipts. A successful face's coefficient and heat requests reuse one resolved
cross query; source resolution allocates and is not an allocation-free kernel
claim. Existing material transport omits cross coupling by default.
Manufactured surface fixtures independently verify the resulting Soret/Dufour
flux and entropy across changing T,p,x, plus source pinning, absent species,
typed admission and atomic refusal; they are not measured-material qualification.

For binary source data, `MaterialBinarySoretDufourTransport<L>` accepts one
`BinaryThermalDiffusionFactorMaterial` with an explicitly ordered species pair,
a finite signed dimensionless property and its original material claim.
The local factor is defined by
`grad ln(x_2/x_1) = alpha_12 grad ln(T)` at zero species flux and constant
pressure (Van-Brunt, Farrell and Monroe, <https://arxiv.org/pdf/2109.05082>,
section 6.3). Species order is never inferred from names or masses. Reversing
the pair requires reversing alpha's sign. An integrated experimental separation
ratio is admissible only when its source supports this local constitutive use;
a Soret coefficient with units 1/K is not silently treated as dimensionless.

At each actual T,p,x query the original source resolver pins and refreshes the
factor claim. The adapter converts it to the same mass-frame kinetic closure:
`H_1-H_2 = alpha_12 R T Mbar/(M_1 M_2)`, using deterministic gauge `H_2=0`.
The original shared face then produces both Soret species flux and reciprocal
Dufour heat, with the original zero-diffusive-flux conductivity convention,
partial-enthalpy accounting and entropy admission. The underlying law still
supplies independent diffusivities, viscosity and conductivity. The adapter
accepts exactly two canonical caloric species including structural zeros; it
does not extrapolate a binary measurement to a multicomponent bank or seed
absent species. Material restrictions govern pure-species endpoint support.
The complete original query, factor and derived heats are available through
`ResolvedBinaryThermalDiffusionFactor`; failed queries leave outputs and the
last complete factor receipt unchanged.

The manufactured public tests in
`tests/material_parcel_transport/binary_factor.rs` check unequal-mass binary
fluxes and conjugate entropy against independent formulas, the measured
zero-flux gradient relation, reversed pair orientation, source interpolation
and claim pinning, structural zeros, unit/assignment/cancellation refusal, and
actual closed-parcel species/energy updates. These are authored regression
checks, not measured-material qualification or a native validation claim.

The shared axial/network/planar/spatial face adds the balanced Soret force
`d_i^T = x_i M_i (H_i-Hbar)/(R T) grad(T)/T`, with `Hbar=sum_i Y_i H_i`,
to the original ideal chemical/pressure forces before calling the same
Maxwell--Stefan owner. Its reciprocal reduced heat flux is
`q_reduced = -lambda grad(T) + sum_i (H_i-Hbar) Ji`.
Here lambda means conductivity at zero diffusive species flux. These kinetic
cross coefficients are distinct from original partial enthalpies: total energy
contains partial-enthalpy diffusion, this reduced heat flux and stress work
exactly once each. A common shift of every H is removed before arithmetic;
zero coupling recovers the original molecular law. Exact absent species are
neither seeded nor assigned spurious flux. At impermeable walls, zero species
flux makes the reciprocal heat contribution zero.

For positive fractions, the local conjugate entropy balance is
`-R sum_i (Ji/M_i) d_chem_i/x_i - q_reduced grad(T)/T²
 = sigma_pair + lambda [grad(T)/T]² >= 0`.
This is the Onsager Schur-complement construction (Van-Brunt et al.,
<https://arxiv.org/pdf/2109.05082>, equations 3.5–3.9). The existing parabolic
engineering scale includes both thermal/composition coupling scales and their
product; no variable-coefficient spectral or fully discrete entropy theorem is
claimed. Bounded scratch remains allocation-free on successful numerical faces.
Focused fixtures independently reconstruct unequal-mass Soret-only, Dufour-only
and combined fluxes and entropy; common-heat gauge/zero-coupling/spatial reduction;
source-axis/domain/cancellation refusal; and periodic composition evolution with
closed original species and total-energy balances. Authored fixture coefficients
are not empirical transport validation.

### Finite material gas inventories with thermal cross diffusion

With `material` and `experimental-euler`,
`material::transport::parcel::ParcelMolecularTransport<L>` consumes the same
prepared `PlanarTransport<L>`, including `MaterialSoretDufourTransport`, on
explicit continuum links between rigid `ClosedGasParcel` inventories. Every
link uses zero barycentric bulk velocity and the original canonical mass-frame
species flux. These are separately exposed continuum paths; they do not add a
second molecular law to the same dusty-gas pore or supply thermal transpiration.
The original caloric species must match exactly, including initially absent
species when they enter a receiving parcel.

Each distinct participating original and proposed endpoint state is qualified
through the shared zero-gradient face kernel, in addition to every actual face
query. The kinetic, material and thermal-cross source domains therefore apply
jointly to the admitted states. The law retains its last successful query,
which may be an endpoint qualification. Per-link reports separately retain the
actual transported species, face enthalpy temperature and reduced heat flux.

All continuum fluxes and native `GasTransportPaths` sample one unchanged
inventory and enter the existing material snapshot reconstruction together.
Each species carries its original donor partial enthalpy once. Added energy is
`q_reduced*A*dt + sum_i dm_i [h_i(T_face)-h_i(T_donor)]`, with the differences
evaluated by the original analytic caloric owner before rounding. Neither large
formation references nor kinetic heat-of-transport datums become extra heat.
The native combined work bound, gross original-species withdrawal, source catalog
join, fixed absolute mass/energy budgets and closed endpoint entropy gate remain
authoritative. Source refusal, cancellation or an overlarge explicit step returns
no candidate state. This is a first-order physical candidate operator; it adds no
separate physical solver. Its whole-window accuracy owner is described below.

`tests/material_parcel_transport.rs` covers independent binary Soret/Dufour
temperature and inventory predictions, unequal-mass frame/orientation,
formation/kinetic-datum invariance, shared nozzle withdrawal, absent species,
initial and final source-domain refusal, and cancellation. Its coefficients are
manufactured integration fixtures, not measured-material qualification.

`ParcelMolecularTransport::{advance_candidate,advance_to}` completes a positive
`[start,target]` interval using bounded full-versus-two-substep comparisons of
those same native physical candidates. Only the fine path enters accepted
states and accounts. Per-parcel species inventory (kg) and temperature (K)
errors use absolute tolerances plus a relative allowance on actual changes,
so large carriers or absolute temperatures cannot conceal a small response.
The fine first-order estimate includes the correction for unequal represented
half intervals. Endpoint clocks define every physical duration; the target is
exact, a shorter final remainder still needs accuracy, and an unresolved
midpoint refuses. An adjacent representable final tail is coalesced into the
proposed span before its physical probes, never by snapping an already solved
clock. The estimate is local and asymptotic, not a global error bound.

The evaluation limit counts every full/fine probe, including discarded and
failed probes. Only typed original withdrawal, pressure-order, entropy or
caloric bracket/root-convergence refusals, and temporal error excess, can reduce
the interval. Original temperatures must lie inside the supplied caloric
brackets. Actual material, kinetic and thermal-cross source/domain refusals are
hard failures. Every physical trial retains initial/final source qualification,
the original enthalpy convention and simultaneous native `GasTransportPaths`.
No extrapolated state, source substitution or clipped inventory is introduced.

Accepted gross species streams and donor partial enthalpies reuse the native
transfer account; additional internal energy retains both directions on each
canonical parcel pair. Rejected probes contribute only to work counts. Detailed
source reports are retained for the final accepted substep. The original material
`audit_closed_window` independently evaluates actual original-to-final total and
source-qualified species mass, NASA internal-energy increments with common
species references cancelled, and caloric/mixing entropy. It and the compensated
sum of accepted native energy residuals must fit the unchanged absolute physical
budgets, with no multiplication by step count. Late refusal or cancellation
publishes neither a gas-map prefix nor an advanced caller clock.

The public adaptive fixtures exercise an independent exact coupled two-mode
Soret/Dufour relaxation, separate species/temperature refinement, stiff pure-gas
Fourier relaxation, accepted-only replay with native nozzle paths, reference
invariance, typed physical retries, late source/work/cancellation rollback and
ordinary plus large represented clocks.

`ParcelMolecularTransport::reservoir_candidate` connects independently exposed
held material baths to finite gas traces through the original zero-barycentric
Maxwell-Stefan/pressure/Soret-Dufour/Fourier face. Every bath reads one original
finite inventory. Signed species mass, combined gross withdrawal and separate
reduced heat enter the existing open species/enthalpy inversion once. Both
stream directions carry original bath partial enthalpy via analytic increments;
there is no fictitious finite bath volume or additional caloric solver.

Original bath, gas, sampled face and final gas source domains qualify. Nonzero
exchange at zero bath species activity refuses without introducing a trace.
Acceptance uses actual endpoint mass, source-qualified species, reference-free
energy residual and gas-plus-bath entropy against fixed absolute budgets. Gas
entropy may decrease when the enclosing bath account compensates it. The
candidate remains first order and requires an admissible caller-supplied step;
it does not add bulk advection, finite bath depletion or an adaptive timeline.
Public fixtures check independent flux/temperature/entropy equations, unequal
molar masses, datum invariance, pure heat, absent species, simultaneous baths,
supported entropy-decreasing overshoot, and late source/cancellation refusal.

`ParcelMolecularTransport::network_candidate` combines internal continuum links
and held-bath faces on that same original snapshot. Its internal face sampler is
shared with the closed finite-gas consumer. All signed species streams, original
donor/midpoint enthalpy corrections and separate bath reduced heat enter one
native gross-withdrawal budget, source join and caloric reconstruction. Bath
arrivals cannot replenish an internal outflow allowance. Internal transfers
cancel in the final open mass/energy/entropy audit; only actual bath streams
cross its accounting boundary. All participating final gas states, including
cells with no direct bath contact, must satisfy the original constitutive law.

Both face families share globally unique physical IDs and one native work bound
before rate allocation. Reports retain internal and bath faces separately;
`reservoir_candidate` delegates with no internal links. Empty bath/internal
families preserve the corresponding physical specialization. This remains a
first-order continuum candidate with an explicit duration. Public network
fixtures cover independent two-cell/bath Fourier balances, simultaneous binary
species/Dufour increments, closed/open parity, shared gross withdrawal, joint
work/identity limits and late source/cancellation atomicity.

`advance_network_candidate` and `advance_network_to` evolve that simultaneous
network across a complete represented time interval. They share the existing
full-step/two-half-step controller with closed molecular transport. The fine
trajectory alone enters accepted states and accounts; no extrapolated state or
rejected probe enters a physical ledger. Local error scales separately bound
each species inventory, temperature, directional species stream and additional
heat stream, so throughflow and reversals cannot hide behind a small net change.
The stream-relative allowance scales transported amounts, not carrier inventory.

Accepted bath streams retain compensated signed sums and separate gross inward
and outward totals. The final original-to-final audit uses unchanged original
bath models, partial enthalpies/entropies and integrated reduced heat; changing
per-step relative caloric charts are not integrated as external energy. Internal
streams cancel at this boundary. Accumulated native caloric defects and actual
endpoint mass/species/energy/entropy each keep one fixed whole-window budget.
Retained accounts share the native work bound, and graph size is admitted before
allocating the held-bath map. Typed physical/time-error refusals refine the step;
source-domain failures and cancellation remain terminal and atomic with the
caller clock. These are local temporal estimates, not certified global errors.
Public fixtures compare two independent thermal modes and a direct bath-heat
integral, isolate stream accuracy from state accuracy, retain balanced-bath gross
species/heat and mixing entropy, and check face-order, work, source, cancellation
and representable-clock behavior. Finite reservoir depletion and native bulk or
dusty-gas path composition remain outside this continuum-plus-held-bath API.


## Reacting Navier-Stokes composition with original-source molecular transport

`ReactingEuler1d::step_with_dissipation` composes exact-duration molecular and
chemistry half steps around conservative advection:
`D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2)`. Each D stage may subcycle its
parabolic update. Only A advances the clock. If heat release changes acoustic
stability and A accepts a shorter duration, all five stages restart from the
accepted field at that duration. Cancellation, source refusal and exhausted
retries cannot publish a partly reacted or diffused field.

The final D endpoint is checked through the original caloric owner against
the reactor's declared temperature interval before publication. Viscous
heating or species enthalpy transport may otherwise leave a narrower chemistry
interval even while the molecular/flow source still admits the state; this
condition triggers bounded complete-step halving without temperature clipping.

The returned optional `SplitDissipationReport` retains both accepted molecular
halves. Whole-step species, momentum, mass and physical total-energy residuals
include both molecular boundary transfers, advection/area terms, and chemical
species/reference-energy changes. No extra chemical or viscous heat source
is added to conservative total energy. The existing `step` method retains
the chemistry/advection/chemistry path with no molecular stages.

Five coupled native tests pass: analytic advection/diffusion/isomerization with
mesh refinement; heat-release-driven duration restarts; closed-wall molecular
traction and zero wall energy work; late cancellation and reused workspaces;
and viscous heating against a narrower chemistry interval, including atomic
refusal and successful shorter replay. The analytic concentration L1 errors
are 6.60993650e-5 / 2.03342713e-5 at 40/80 cells. The complete isolated native
source-linked suites pass 132 thermochemistry and 79 gas unit tests; these
results are scoped to those owner builds, not the full runtime dependency graph.

`cargo run -p fs-gas --features experimental-euler --example reacting_oxygen -- 40 molecular`
enables the complete original-source path in the O/O2 example. Calorics,
reaction and Lennard-Jones parameters retain Cantera v3.2.0 GRI3 attribution
and original blob 9a62004ec17bbec8ac11e97785932d426dd58610. At 40 cells over
1 ms the native run accepted 113 coupled steps / 226 molecular substeps,
with mass residual -4.336809e-19 kg, physical total-energy residual
-6.274981e-13 J, and zero reported momentum residual. Atomic oxygen mass
falls from 4.237578e-5 to 3.380238e-5 kg. These are numerical consistency
checks of an isolated reaction, not an experimentally validated air mechanism.

The transport law's dilute-neutral-gas, Eucken-conductivity, explicit bulk
viscosity and constant-area restrictions remain in force. This symmetric
composition and either selected chemistry method do not by themselves
establish global second-order accuracy or a fully discrete entropy inequality.


### RKL2 diffusion in reacting flow and complete physical windows

`ReactingEuler1d::step_with_stabilized_dissipation` opts both molecular halves
of the original D-C-A-C-D composition into the shared RKL2 owner. Chemistry,
advection, caloric support, final chemistry-domain admission, complete-step
restarts and physical boundary/source balances remain with their original
owners. The default explicit APIs retain their original composition.

A single actual RHS counter spans both halves and every discarded coupled
attempt. Exhaustion is terminal; halving the physical duration cannot create
another allowance. Reports retain the original coupled balance and both
accepted RKL2 half reports, while the aggregate RHS count includes rejected
work.

`advance_stabilized_to` completes an explicit represented target through that
same step owner; `advance_stabilized_for` admits a duration only when adding it
to the clock preserves its span within a named 1e-10 relative roundoff bound.
Every terminal remainder is integrated. Actual accepted durations and physical
source/boundary transfers accumulate with compensation, and integrated time
must agree with the represented window. Coupled attempts and molecular RHS
calls have separate whole-window limits. One lazily reused conserved-state
backup restores the original accepted field and clock after any late refusal,
cancellation, budget exhaustion or final balance failure.

This extension retains stationary constant-area molecular geometry and
original source-owned chemistry/transport closures. Macrostep and chemistry
controls still govern accuracy; the interval owner supplies no global temporal
estimator, nonlinear RKL2 stability certificate, moving-mesh coupling or
experimental material validation.


Seven new focused regression groups pass with the original reacting and
molecular tests: 16 reacting tests and 33 dissipation/planar/spatial tests,
49 total. The native selected-owner build uses the genuine current
Maxwell-Stefan v3 thermochemistry and original gas/chemistry physical modules;
this is not a full workspace Cargo gate. The independent advected
reaction-diffusion mode has L1 errors 2.61263736e-5 and 7.60930020e-6 at
40 and 80 cells. In the authored 96-cell comparison, a 2.5e-4 s molecular
macrostep cap uses 36 actual RHS calls versus 132 explicit SSP RK2 calls,
with maximum/mean conservative-component discrepancies 7.29747e-6 /
2.36738e-7. These are fixture-specific operation counts and numerical
accuracy checks, not measured runtime or material-validation claims.

### Stationary reservoir baths through reacting fields

`ReactingEuler1d::step_with_reservoirs` advances the existing C-A-C composition
against prescribed species-bearing exterior cells.
`step_with_dissipation_and_reservoirs` and
`step_with_stabilized_dissipation_and_reservoirs` pass those SAME bath states to
the advective stage and both molecular halves of the original D-C-A-C-D owner.
Every whole-step retry starts from the accepted field and replays the baths
over the actual shortened duration. Chemistry advances only the interior.
Original characteristic inflow/outflow, half-cell molecular gradients, bath
enthalpy and stationary traction work retain their constitutive owners.

`advance_stabilized_with_reservoirs_for` and
`advance_stabilized_with_reservoirs_to` extend the same coupling to complete
physical windows. Original compensated advective, molecular, chemical and
left/right boundary accounts include each transfer once. A single molecular
RHS budget covers both halves, all retries and all accepted steps; the original
separate coupled-attempt bound and complete-window rollback remain in force.
Baths are admitted even on a zero-duration window. Two absent baths preserve
the original state, reports and cancellation sequence.

Bath cells must use the original canonical species order and energy gauges;
their borrowed representation has no source identity token. Prescribed baths
are stationary and undepleted, and periodic fields refuse them. The existing
finite network/plenum owner supplies finite storage. Moving baths, general
boundary histories within one call and a global splitting-error estimator are
outside this interface. The new focused regressions cover an independent
Dirichlet reaction-diffusion eigenmode, characteristic inlet species/enthalpy,
heat-release retries and shared work budgets, full-window rollback/reuse and
absent-bath equivalence. These are numerical/source consistency tests, not
experimental reaction or transport validation.

## Characteristic ports with distinct local compositions

`MixtureGas::reservoir_flux` connects a conservative duct trace to a stationary
species-bearing reservoir in the same original source bank and per-species
energy gauges. Local compositions may differ. One shared borrowed-caloric
solver now supplies both frozen and mixture ports: the duct composition owns
the outgoing pressure wave, while the reservoir composition owns isentropic
acceleration and the sonic inlet limit. The advected donor supplies every
species flux and its total enthalpy; structural species zeros remain zero.

The public output buffer is updated only after the complete port has passed
source, numerical and cancellation checks. Successful evaluation uses bounded
stack scratch and no replacement gas owner or heap workspace. A strong inflow
whose auxiliary velocity bound reaches the full cooling interval locates its
supported sonic throat before constructing an unused cold endpoint; actual
port admission remains strict.

Five independent regressions cover two different gamma laws across a port,
analytic choked inflow and sonic outflow, branch acoustic impedance, donor
enthalpy, original formation-reference/orientation invariance and unchanged
output on late cancellation. All 72 focused staged gas tests passed, including
the existing frozen networks, wall reflections and reacting-flow checks.
This run precedes the concurrently added molecular-dissipation module and is
not a full workspace Cargo/DSR verdict. The closure remains an inviscid compact
reservoir model, without resolved junction mixing or empirical discharge loss.

## Conservative networks with evolving species

`multicomponent::network::MixtureNetwork1d` connects the existing conservative
mixture branches to finite rigid `Plenum` volumes. Every owner shares the exact
original species bank, canonical order, caloric path and individual energy
references; its local composition and temperature may differ. The shared
characteristic port supplies donor species, pressure impulse and total
enthalpy once. The corresponding branch-end and plenum accepted exchanges
are exact opposites for each species and energy.

One global SSPRK2 transaction advances all branches and plenums. Actual port
speeds constrain branch and compact-volume CFL at both stages; every retry
restarts all owners. State, clocks and reports publish together after final
admission/cancellation. Flat numerical/report storage is prepared once and
reused, and compact-volume temperature inversions use cached guesses. Compact
storage has no bulk momentum; its fixed support receives the reported local
port impulse while incoming kinetic energy remains in the energy balance.

Four focused regressions cover different-composition rest in tapered geometry,
finite-pulse splitting and actual species transport, independent extensive
species/physical-energy accounting, accepted port exchange, compact-volume CFL,
source/topology refusal, cancellation/reuse, orientation and formation gauges.
The integrated native gas build, including the concurrently landed molecular
and reacting transport modules, passes 88 tests.

`cargo run -p fs-gas --features experimental-euler --example mixture_manifold -- 64`
runs pressurized N2 into an air surrogate containing N2/O2 through three passages
and a 30 mL plenum. It retains the original 200--1000 K NASA-9 coefficients from
Cantera v3.2.0 `data/airNASA9.yaml`, blob
`5975d1303819332c1bd77a783d9873b66f65c952`, with the declared 1 atm reference.
The 16/64-cell runs reach 4 ms in 149/593 steps, with no HLL fallbacks or
halvings. Physical-energy residuals are -3.170523e-14 / 1.919073e-15 J;
the largest species mass residual is 1.084202e-19 / 1.355253e-20 kg.
Both daughters receive nitrogen and the final plenum N2 mass fractions are
0.973409707 / 0.979168255.

All sampled states remain inside physical species bounds and the source path.
A small composition-extremum undershoot remains in this conservative MUSCL
baseline (64-cell daughter minimum N2 fraction about 0.766145 versus initial
air 0.767082); no composition maximum principle or mesh-converged prediction
is claimed. This network transport owner has no chemical source, resolved
junction mixing, molecular port diffusion, wall loss or entropy certificate.
The required `dsr quality --tool frankensim` could not run because DSR is
unavailable; focused native checks are not the workspace gate.

## Finite-rate reacting shockwave networks

`multicomponent::network::reactive::ReactingNetwork1d` owns a complete mixture
network and one prepared original-source isochoric reactor. Every branch cell
and finite rigid plenum participates in `C(dt/2), A(dt), C(dt/2)`, using its
actual volume. The existing reacting single-branch owner and this network
share the same per-cell chemistry update. Momentum and physical chemical total
energy are fixed during C; the stored referenced-energy change is exactly the
negative change in species reference energy, without a second heat-release
source. Endpoint sensible energy is independently reconstructed from the
original species models.

The initial duration respects branch and compact-volume acoustic scales. If
the first chemical half raises wave speeds and transport accepts less time,
the whole composition restarts from accepted branch and plenum states at the
actual duration. Exhausted retries, source refusal and cancellation cannot
publish a partly reacted network. Fields, clocks, combined reports and port
exchanges publish together after final admission. Trial fields, amount storage,
reactor and numerical workspaces are reused; reset work polls in bounded tiles.

Combined reports distinguish chemical species/reference-energy changes from
external transport and reconstruct physical total-energy balance. Accepted
port exchanges remain opposite to the branch transport reports. Those reports
explicitly describe the intermediate transport stage; releasing the raw network
clears its stale transport-only whole-step report. Compact plenums retain the
same perfectly mixed, zero-bulk-momentum, fixed-support interpretation.

Four new regressions compare homogeneous reversible isomerization with analytic
decay, and strong exothermic steps with independent constant-volume reactors at
the accepted duration after a heat-driven restart. They separately reconstruct
branch/plenum chemistry and actual port-energy budgets, reject mismatched
original sources and plenum temperatures outside the chemistry interval, and
verify late cancellation plus workspace reuse. The complete integrated native
gas suite passes 92 tests, including existing reacting molecular-transport
checks. Eight existing nonlinear duct tests pass against the shared port and
mixture-network kernel. These are production-source native checks, not a full
workspace Cargo/DSR result; DSR is unavailable in this environment.

This owner composes inviscid network transport with finite-rate chemistry.
Molecular junction diffusion, wall heat/loss, valves and resolved mixing remain
separate closures. The chemistry entropy diagnostic is not a network entropy
certificate, and the selected chemistry method combined with symmetric
splitting does not establish global second-order convergence.


## Two-layer compartments with conservative moving-interface work

`multicomponent::compartment::TwoZone` adds a reusable two-zone control volume
for issue #2. It retains the original `MixtureGas` species bank and source-path
admission. Each layer stores canonical species masses, sensible internal energy
in that bank's gauge, physical volume and temperature. A layer may initially
be absent exactly; real entrainment and heat can create it without artificial
mass, species or temperature seeds.

A source increment carries heat and advected partial enthalpies. The fixed room
volume is partitioned by common pressure and the actual variable-Cp ideal-gas
EOS. The implicit midpoint interface work is
`W = (p_old+p_new)/2 * (V_lower_new-V_lower_old)`; the other layer receives
exactly `-W`. A bounded scalar energy-partition solve recovers both temperatures
through the original gas inversion. Its derivative is positive on the admitted
positive-Cv source paths. Neither fixed gamma nor a new caloric surrogate is
introduced. The total sensible-energy source is conserved by construction;
physical-energy reports retain the original species formation references.

`inflow` and `carried_rate` provide canonical enthalpy-bearing source rates.
The same transfer is added to the receiver and negated for the donor.
`advance` completes exactly one supplied-source duration atomically; a network
or time-quadrature owner must recompute its rates when retrying. Prepared rooms
can copy accepted physical state only when their source and geometry agree.
Cancellation, negative inventory, unsupported calorics, an out-of-room interface,
or a bounded solve failure publishes nothing and leaves the clock unchanged.

`TwoZone::apply_split` composes a supplied heat/stream source substep with an
already timed outer update. It uses the same atomic inventory, caloric,
pressure and interface-work solve as `advance`, then restores the exact
original clock after success. Only the enclosing time owner advances the
physical interval and records those source rates once. A focused regression
compares its complete physical state and formation-energy-aware stream balance
with an ordinary source advance while checking the retained clock and atomic
cancellation. This operation does not add chemistry or a higher-order split.

Seven native tests pass: the constant-Cp CFAST pressure reduction and adiabatic
interface refinement; independently integrated variable-Cp heating; conservative
unequal-species enthalpy exchange; formation of an initially absent upper layer;
original EOS, composition and interface checks; cancellation/refusal/reuse;
and representability admission including nonbinary room volumes. The governing
mass, energy, ideal-gas and fixed-volume equations are documented in
[NIST TN1889v1, chapter2](https://doi.org/10.6028/NIST.TN.1889v1).
This is an Estimated well-mixed low-Mach compartment owner, not a flame,
ceiling-jet, acoustic, soot-radiation or experimental-validation claim. Prescribed
source rates are frozen over a step; the midpoint work alone does not make the
complete source/network integration second order.


## Inventory-limited compartment burning

`multicomponent::compartment::burning` binds an original admitted ideal-gas
kinetics mechanism and one conserved global reaction to the exact room species,
NASA calorics, and per-species energy references. The requested nonnegative
forward extent is capped by every net reactant inventory, including oxygen,
with optional nonnegative incoming species rates over the explicit interval.
Incoming rates are availability only: the caller must add their physical mass
and enthalpy exactly once. A driver with simultaneous reactant outflow must
reserve that withdrawal first or transport before burning; timestep halving
alone cannot make an additional withdrawal compatible with complete depletion.

The adapter returns chemistry-only species rates and the source-coefficient
energy-reference correction `-extent_rate * Delta u(Tmin)`, evaluated at the
original `MixtureGas` reference temperature. Conserved elemental formation
references cancel before absolute-energy rounding; distinct formation energies,
original caloric corrections and molecular-count pressure work remain active.
Reported reaction enthalpy and power are diagnostics, never a second heat
source. Power-to-extent conversion uses source-coefficient Delta h at an explicit
reference temperature. `chemical_reference_roundoff_w` separately exposes the
discrepancy from independently rounded species-reference powers and is never
added as heat. The original moving-interface solve retains all paired p dV work. One common extent may be rounded down by a bounded
number of adjacent floating-point values to admit actual update arithmetic;
individual species are never clipped, seeded, or independently repaired.

This is an estimated global, inventory-limited burning closure, not an
Arrhenius ignition/extinction model, equilibrium solver, soot model, or proof of
second-order temporal accuracy. Reactions with species on both sides are
explicitly refused by this single-pass adapter. Six focused tests cover oxygen
exhaustion, atom/mass conservation, analytic variable-Cv formation-energy
heating, inflow counted exactly once, a one-ULP depletion case, temperature-
dependent reaction enthalpy, cancellation, and source identity/refusals.


### Clock-preserving transport/burning composition

`InventoryLimitedReaction::apply_split` evaluates the actual transported
inventory with no repeated incoming source, delegates chemical species/energy
and pressure-interface relaxation to the original atomic TwoZone advance, and
then restores the physical clock infallibly. The source duration controls the
reaction increment but does not advance elapsed time twice. This first-order
split allows prior ventilation withdrawals to reduce oxygen availability before
burning. The returned chemical source is recorded once by a composing network;
reported reaction heat is still never an extra energy source. All eight
original burning tests pass, including two added regressions for post-outflow
oxygen depletion, fixed clock, paired interface work, physical-energy balance,
late cancellation, reuse, and caloric-domain refusal of excessive source heating.

## Molecular transport through finite-plenum networks

`multicomponent::network::dissipation::NetworkDissipationOperator` extends the
existing original-source Maxwell--Stefan, Fourier and longitudinal-viscous
operator across branch/plenum connections. Every branch retains the existing
stationary constant-area restriction. A caller supplies one strictly positive
plenum-side `Length` per connection. Together with the branch half-cell width
it defines the constitutive gradient and interpolated face state; no mixing
length is inferred from plenum volume.

Each connected face uses the same molecular kernel as an interior branch face.
Species diffusion is mass-averaged, carries each original partial enthalpy,
and shares total-energy flux with conduction and viscous work. The branch and
plenum receive equal-and-opposite species and energy exchange. A stationary
plenum retains zero bulk momentum; its fixed support receives the signed local
axial impulse. No scalar momentum is summed over differently oriented branches,
and no separate viscous-heating source duplicates conservative energy.

One global explicit SSPRK2 subcycle advances all participating physical
volumes. Its estimated parabolic bound includes both branch cells and the sum
of all port rates divided by each plenum volume. Every requested molecular
interval completes its exact duration before fields publish; the flow clock
is unchanged. Law, constitutive solver, stages and diagnostic arrays are prepared
once and reused. Candidate port/branch reports remain separate from accepted
reports, allowing a composing owner to reject an entire split step atomically.

`multicomponent::network::dissipating::DissipatingNetwork1d` owns a complete
`D(dt/2), A(dt), D(dt/2)` step. If advection shortens the duration, both molecular
halves restart from the original accepted state. Its report closes species and
physical energy against the combined external boundary transfers, and its
accepted molecular-port slice is the sum of both D halves. The A-only exchange
and branch reports explicitly describe the intermediate advection stage.
Releasing the raw network clears its stale intermediate whole-step report.

Five molecular regressions check two-volume binary/Fourier exponential
relaxation and timestep refinement, source-consistent partial-enthalpy transfer,
viscous support impulses, formation-gauge/orientation invariance, compact-volume
subcycling and atomic cancellation/reuse. Three composing-owner regressions
check zero-law identity, independent A+D port/storage budgets, late cancellation
and exhausted molecular work without partial publication. All eight new checks
passed in the integrated production-source native build.

The compact gradient-length closure is an Estimated continuum model. It does
not resolve multidimensional mixing, radial wall layers, variable-area viscous
stresses or rarefied gas. Thermal cross coupling requires the explicit source
law above; a fully discrete entropy or invariant-domain theorem remains outside
this operator. Native source-linked
checks are supplementary: `dsr quality --tool frankensim` is unavailable here
(`dsr: command not found`), so no workspace DSR verdict is claimed.

### Source-conjugate molecular entropy through finite networks

`NetworkDissipationOperator::advance_entropy_consistent_for` applies the original
axial source-conjugate molecular kernels at both branch interiors and finite
plenum ports. Divided source enthalpies and logarithmic chemical efforts retain
the original Maxwell--Stefan, Fourier, reciprocal Soret/Dufour and viscous-work
laws. Actual plenum mixing lengths define gradients. Both finite storage owners
receive opposite species/energy transfers; fixed plenum supports receive axial
impulse without mechanical work or a second viscous-heating source.

One original-source entropy account covers all branch cells and finite plenums.
Every synchronized SSPRK2 substep and the direct original-to-final interval
must satisfy the same declared entropy allowance and individual species
continuity bounds before publication. Internal ports supply no exterior entropy.
Only the source-reference contribution of numerical species drift is excluded
from admission; full source entropy is reported separately. Native candidate
temperature inversions are reused, with no flattened original field copy or
per-substep allocation. The flow clock remains unchanged. Final refusal,
exhausted work and cancellation preserve all accepted fields and receipts.

`advance_for_with_entropy` measures the original engineering faces with this
same account and retains their accepted native trajectory. It can refuse an
energy-conserving face with a negative spatial entropy rate. The conjugate entry
requires each species to be positive on both sides of a face or absent on both;
one-sided absence refuses without a floor or alternate flux. These closed,
stationary constant-area entries exclude external baths, moving volumes,
chemistry, advection and RKL2 composition. Measured finite-step admission is not
a fully discrete entropy or invariant-domain theorem. Focused G0/G1/G3/G4
regressions exercise negative-entropy repair, variable-Cp storage and support
work, orientation/reference invariance, native-trajectory parity, the complete
interval allowance, species support and final-poll cancellation/reuse.

### Stabilized molecular intervals with finite plenum storage

`NetworkDissipationOperator::advance_stabilized_for` opts the original complete
network operator into second-order Runge–Kutta–Legendre (RKL2) macrosteps.
`advance_stabilized_with_reservoirs_for` retains prescribed exterior baths on
the same path. Both reuse the axial Legendre coefficients and branch recurrence,
the existing network spatial kernel, and the original source/geometry admission
and final publication owner. Each finite plenum evolves with its actual volume
and the sum of every incident port flux. Branch fields, plenum inventories and
signed endpoint transfers follow the same recurrence, including its negative
weights, so internal species and energy exchange remains equal and opposite.
The stationary plenum's support retains each port's local axial impulse.

`StabilizedDissipationControl` retains the original rate and residual controls.
The global rate includes all branch cells and every finite-plenum port rate
divided by its volume. For `s` stages the estimated forward-step interval is
extended by `(s*s+s-2)/4`; the original source, caloric, constitutive and global
rate checks run at every intermediate state and the final endpoint. Refused
trials halve the whole macrostep and restart from its cached original RHS.
Stage, attempt, halving and cumulative actual-RHS budgets are bounded; an
accepted `s`-stage macrostep uses `s+1` spatial evaluations. The report separates
the unscaled parabolic number, its stability-extension-scaled value, the largest
accepted stage count and actual work, including discarded attempts. Composing
owners can carry that same work counter across both molecular halves and
complete-step retries without publishing trial reports.

The three existing conservative state arrays are reused. Anchor-face,
anchor-plenum-RHS and endpoint-recurrence buffers are reserved lazily and reused
without per-stage allocation. Admission, signed branch/port/exterior balances,
and complete-duration publication are shared with SSPRK2. Cancellation, budget
exhaustion or source refusal preserves all accepted fields, clocks and reports.

The physical macrostep cap remains an explicit accuracy choice; there is no
temporal error estimator, unconditional nonlinear positivity guarantee or
entropy certificate. Authored regressions independently evaluate a finite-volume
two-pool Legendre eigenmode with a stiff small plenum, Fourier relaxation and
temporal refinement with unequal heat capacities, and an affine three-storage
bath/plenum matrix polynomial. They also check actual stability work versus SSP,
signed endpoint and finite-inventory transfers, zero net diffusive mass,
orientation and formation-gauge invariance, late cancellation, cumulative work
refusal and reuse. The work comparison is a stability comparison, not a claim
of equal temporal error. Test execution is not implied by this contract.

### Source entropy for stabilized finite-network intervals

`NetworkDissipationOperator::advance_stabilized_entropy_consistent_for` connects
the original source-conjugate branch/finite-plenum faces to the original RKL2
recurrence. Divided caloric enthalpies, logarithmic chemical efforts,
Maxwell--Stefan diffusion, Fourier heat, authored reciprocal cross effects and
viscous work remain in their existing owners. All finite inventories and signed
port transfers retain the same Legendre weights. Internal ports supply no
exterior entropy and stationary plenum support performs no mechanical work.

Every physical macrostep endpoint and the direct complete interval must pass
the original source-entropy and individual-species continuity bounds. Internal
Legendre states retain native source, caloric, constitutive and estimated-rate
checks; they are not independent physical entropy-producing substeps. Negative
endpoint entropy invokes the original bounded halving retry from its unchanged
anchor. Native endpoint temperatures are reused, and branch/plenum volumes are
borrowed directly without another EOS inverse, flattened storage copy,
per-stage allocation or extra constitutive RHS.

This entry admits stationary closed constant-area geometry with common species
support at active molecular faces. Baths, moving volumes and other sources are
separate paths. The physical macrostep cap remains an explicit accuracy choice.
Linear RKL stability does not establish nonlinear positivity, a fully discrete
entropy theorem or global temporal order. Focused G1/G4 regressions compare a
stiff finite-plenum Fourier mode with an independent Legendre polynomial and
source-entropy integral; compare actual stability work with SSP without claiming
equal temporal error; retain unequal calorics, common absent species and
formation-reference invariance; and refuse a negative-entropy RKL endpoint
despite positive smooth conductivity and admitted intermediate states.

### Accuracy-controlled molecular intervals with finite storage

`NetworkDissipationOperator::advance_adaptive_for` and
`advance_adaptive_with_reservoirs_for` add explicit physical local-accuracy
controls to this same network RKL2 owner. `AdaptiveNetworkDissipationControl`
retains the original source, stage, attempt and complete-window RHS limits and
adds absolute species-density, temperature and axial-velocity accuracy plus a
componentwise relative tolerance. Each trial compares one full step with two
half-steps using the SAME selected stage count, hence the same second-order
method coefficients. The difference divided by three estimates the local error
of the two-half result. Selecting fewer stages for the halves would change the
leading error coefficient and is deliberately not done.

The maximum error norm checks every original species separately in every
branch cell and finite plenum, together with temperature and velocity decoded
by the original caloric owner. Large total mixture mass or arbitrary formation
references cannot hide trace-species or thermal error. Every individual signed
endpoint transfer is also checked. Its absolute species accuracy uses the
smaller incident storage volume, thermal accuracy uses the smaller actual
incident heat capacity, and support-impulse accuracy uses the branch-cell mass.
Relative transfer accuracy scales the actual transferred quantity; cancellation
between two ports or between branches does not erase a quadrature error.
Canonical sensible energy is used for this port norm, retaining formation-gauge
invariance. All actual physical energy receipts retain their original meaning.

An accepted interval publishes the admitted two-half field and its matching
conservative endpoint quadrature, without Richardson extrapolation, clipping or
post-hoc balance repair. Failed accuracy probes shrink with a bounded
second-order local-error controller; failed physical probes retain the original
halving policy. Two field arrays and two endpoint arrays are reserved lazily
and reused, with no per-probe numerical allocation. The initial RHS is shared
by the full and first-half probes, and the admitted midpoint RHS starts the
second half: an unrefused `s`-stage probe costs `3*s+1` spatial evaluations.
All refused probes consume that same whole-window work budget.

`AdaptiveNetworkDissipationReport` nests the original stabilized/balance report,
the largest ACCEPTED component/port error ratio (at most one), and the number of
accuracy refusals. Accepted steps count paired half-step intervals; rejected
steps include physical and accuracy refusals. A composing owner can retain one
RHS counter across molecular halves, geometry changes between those halves,
and discarded complete coupled attempts. Final publication and all accepted
report preservation use the original shared transaction.

This is an estimated local-error lane, not a certified global, spatial or
operator-splitting error bound, nor an unconditional positivity/entropy theorem.
Focused regressions cover a stiff finite-plenum two-half Legendre mode,
tolerance refinement against continuous species and unequal-capacity Fourier
exponentials, affine three-storage bath evolution and port quadrature against a
matrix-exponential reference, formation-gauge/orientation invariance, and late
work/cancellation rollback with reuse. Test execution is reported separately.

### Adaptive source-conjugate molecular entropy

`NetworkDissipationOperator::advance_adaptive_entropy_consistent_for` combines
the original full/two-half RKL2 accuracy controller with source-conjugate faces
throughout every branch and finite-plenum port. Species density, temperature,
velocity and signed port-transfer errors retain their existing componentwise
norms, original caloric inversions and fixed-stage comparison. The accepted
two-half state and conservative transfers must also pass the original entropy
and individual-species continuity gate from the complete probe's unchanged
anchor. The discarded coarse estimate and internal midpoint/stage states do not
constitute separately admitted physical entropy intervals.

Temporal refusals use the native error controller; entropy refusals halve the
complete probe. All discarded work consumes the same RHS budget, while native
endpoint temperatures and the retained adaptive anchor supply entropy without
another EOS inverse, constitutive RHS or per-probe allocation. The direct
original-to-final entropy check precedes the shared atomic publication, so late
refusal preserves accepted fields, clocks and receipts. Reports retain both
local-accuracy diagnostics and measured source entropy.

This entry admits stationary closed constant-area networks and common species
support at each active face. Baths, moving geometry, reactions and advection
remain separate. Its estimated local molecular error and measured finite-step
entropy do not certify global temporal error or an invariant domain. Focused
regressions cover Fourier tolerance refinement, repair of conservative
negative-entropy engineering mixing, unequal molecular masses and variable Cp,
finite-port balances, formation-reference invariance and late rollback/reuse.

### Exterior molecular baths through complete network transactions

`NetworkDissipationOperator::advance_with_reservoirs_for` supplies one optional
left/right prescribed bath pair per branch. Only unconnected physical ends may
use these states; prescribing a bath on a finite-plenum connection is refused.
The original axial admission and reservoir flux kernel retain source calorics,
physical end half-widths, bath species enthalpy and zero stationary-bath
traction work. Bath primitives and volumetric heat capacities are admitted
once per interval into construction-time scratch; original cells remain
borrowed, with no new per-step numerical allocation. As with the characteristic
reservoir API, raw `Cell` values require caller-consistent canonical species
order and energy gauges; they carry no source identity token.

Both `DissipatingNetwork1d::step_with_reservoirs[_and_sources]` and
`ReactingNetwork1d::step_with_dissipation_reservoirs_and_sources` pass the same
baths through both molecular halves and characteristic advection, including
every coupled retry. Previously those prescribed states reached only advection,
leaving the end molecular condition unchanged. The complete transaction now
admits diffusion driven by an exterior composition or temperature difference
even at zero pressure-driven flow. Advected donor enthalpy and molecular bath
enthalpy remain separately accounted and enter the complete balance once each.
Finite-plenum molecular storage and equal-and-opposite internal port exchanges
continue through their existing owner and declared mixing lengths.

Invalid topology/source states, work exhaustion and late cancellation leave
all accepted fields, clocks, branch reports, port exchanges and whole-step
reports unchanged. All absent baths retain the original molecular path and
cancellation sequence. Six focused regressions cover an independently formed
three-storage SSP update with a prescribed exterior bath and finite plenum,
orientation and formation-gauge invariance, Fourier heat/storage balance,
isobaric D/A/D diffusion and complete rollback, bath replenishment of a species
consumed by reversible chemistry, and absent-bath/refusal/reuse behavior. They
are authored physical tests; native execution is not implied by this contract.

### Stabilized molecular transport in complete network steps

`DissipatingNetwork1d::step_stabilized[_with_sources]` and
`ReactingNetwork1d::step_with_stabilized_dissipation[_and_sources]` select the
same prepared network RKL2 operator for both molecular halves. One internal
dispatch serves the original D/A/D and D/C/A/C/D transactions, preserving
characteristic advection, chemistry, source ordering and physical accounting.
The nonreacting controls explicitly select RKL2; no unused explicit molecular
control silently overrides its macrostep or work limit.

Each complete step borrows one molecular RHS counter across both halves and
all discarded acoustic/heat-release retries. Exhausting that allowance is a
terminal original `OperatorWorkLimit`, never a fresh retry budget. Successful
reports retain both accepted RKL2 reports and the complete actual RHS count,
including work from rejected compositions. Original accepted molecular port
receipts sum both halves; exterior diffusion, advected donor enthalpy, wall heat
and prescribed chamber pressure work remain distinct physical accounts.

The source-aware entries admit optional prescribed baths and finite wall/valve
sources. The reacting entry also admits original prescribed plenum volume
rates, so each molecular half evaluates its actual storage volume. The original
moving-advection owner supplies mechanical pressure work once. Any failure
preserves accepted gas, clocks, geometry, finite walls and separately owned
source/operator reports; only private numerical work may change.

Five focused regression groups cover the independent three-storage bath ODE,
finite-plenum uptake and complete balances, both-half RHS exhaustion and late
cancellation/replay, continuous two-pool reaction/diffusion decay, heat-release
retries with one work allowance, and moving chamber/finite-wall energy closure.
These are authored numerical tests, not measured-material qualification.
RKL2's extended estimated stability interval does not certify the complete
composition's time accuracy, positivity or entropy; caller-selected maximum
macrosteps, original stage admission and refinement remain necessary.

### Adaptive molecular transport in complete network steps

`DissipatingNetwork1d::step_adaptive[_with_sources]` and
`ReactingNetwork1d::step_with_adaptive_dissipation[_and_sources]` select the
original network RKL2 full/two-half error controller inside the existing D/A/D
or D/C/A/C/D transaction. The same physical species-density, temperature,
velocity and conservative port-transfer tolerances apply to both molecular
halves. Each half returns its accepted local error estimate and temporal
rejections alongside the original physical balance. No separate transport,
chemistry or source implementation is introduced.

One RHS counter spans full and half probes, both accepted molecular intervals,
failed temporal probes and every discarded complete coupling attempt. Work-limit
exhaustion is terminal; coupling retries cannot reset the allowance. Only the
accepted two-half molecular solutions and their original integrated exchanges
enter branch/plenum inventories and the complete-step species/energy balance.
Temporary estimator paths cannot publish port transfers or source receipts.

Prescribed baths still reach molecular transport and characteristic advection.
Moving-plenum molecular intervals use their actual pre/post-advection volumes;
the existing moving-volume owner applies pressure work once. Finite walls,
valves, gas fields, geometry and accepted diagnostics retain the original atomic
publication boundary. The molecular estimator does NOT control or certify
chemistry, source, advection, spatial, splitting or global time error, and does
not infer entropy stability for their composition.

Focused G1/G4 regressions compare bath transport with an independent continuous
three-capacity ODE, resolve a commuting reaction/diffusion decay without a
manually tuned molecular step, account for discarded heat-release retry work,
close moving gas-plus-wall energy, and preserve accepted fields/geometry/reports
under a second-half work refusal and late cancellation with replay.

## Reacting networks with molecular transport

`ReactingNetwork1d::step_with_dissipation` composes
`D(dt/2), C(dt/2), A(dt), C(dt/2), D(dt/2)` for every branch cell and finite
plenum. It reuses the same prepared molecular network operator and original
constant-volume chemistry update; each molecular half may subcycle, while
only A advances time. If heat release or diffusion changes the admissible
acoustic duration, all five stages restart from accepted physical storage at
the actual duration.

Molecular operator reports remain unpublished through coupled retries and
cancellation. Complete fields, clocks, reaction/transport/molecular reports and
the sum of both molecular port transfers publish only after the final check.
The original C/A/C method remains available without molecular stages. Combined
balances include both external molecular transfers and original chemical
species/reference-energy changes; no duplicate heat-release or viscous-heating
source is introduced. Final molecular temperatures in every branch cell and
plenum must remain inside the mechanism's possibly narrower caloric interval.

Three new coupled regressions pass: a two-pool analytic reversible
reaction-plus-diffusion exponential, heat-release-driven complete-duration
restart with homogeneous chemistry identity, and late cancellation preserving
both network and separately owned molecular reports followed by successful
reuse. All 110 staged native gas tests passed, including the earlier
single-branch molecular/reacting paths and the original compartment tests,
before the subsequent Maxwell--Stefan trace-gradient scaling correction.
After that correction, all ten exact-production-module Maxwell--Stefan tests
and rebuilt 16/64-cell molecular manifold consumers passed. This scoped
validation does not claim tests added concurrently elsewhere on main or the
unavailable workspace DSR gate.

The selected chemistry method and symmetric composition do not by themselves
establish global second-order accuracy or a discrete entropy certificate. The
molecular operator's constant-area, continuum, explicit mixing-length and
transport-law applicability limits remain unchanged.

### Complete source entropy in closed reacting networks

`ReactingNetwork1d::step_with_entropy` admits C(dt/2)-A(dt)-C(dt/2) through the
original complete-network retry transaction. `step_with_dissipation_and_entropy`
adds source-conjugate molecular halves D(dt/2) on both sides. Every branch cell
and finite plenum retains the original fixed-U,V reactor, caloric data, reaction
mechanism and energy-reference conversion. Advection uses the existing measured
network entropy gate; it remains the native engineering flux. Molecular faces
use the existing divided caloric/chemical pairing and authored port lengths.

One immutable pre-step source basis admits the actual original-to-final entropy.
The account observes actual represented chemical species changes per cell and
per half, retaining their physical source-entropy constants. Only numerical
species-source continuity defects lose their arbitrary reference contribution;
each defect is independently bounded against original mass plus absolute
chemical throughput. Internal plenum transfers are enclosed in the state
account. Summing component entropy reports or reactor tolerances cannot replace
the complete endpoint gate, and no reaction or viscous heat is added twice.

Acoustic shortening or entropy refusal restarts the complete composition. Both
molecular intervals finish their exact half duration, but their accepted port
and branch receipts remain unpublished until the complete source/energy/entropy
checks and final cancellation poll succeed. The final molecular state is also
admitted against the original chemistry temperature interval. Refusal preserves
all accepted network fields, clocks and participating operator reports.

These entries admit only fixed closed networks, with stationary reflecting
exterior ends or isolated periodic branches. The molecular entry additionally
requires constant-area branches and common species support at active faces.
Baths, moving volumes, wall/valve sources and stabilized/adaptive molecular
composition are excluded. Measured component and whole-step admission does not
establish a fully discrete coupled entropy or global temporal-order theorem.
Focused G1/G4 regressions reconstruct source entropy with variable Cp and a
nonzero physical species entropy datum, create an initially absent product,
exercise heat-release duration retries, and cancel/reuse the final transaction.

### Complete reacting entropy with stabilized molecular transport

`ReactingNetwork1d::step_with_stabilized_dissipation_and_entropy` uses the
source-conjugate RKL2 molecular halves inside the original D/2-C/2-A-C/2-D/2
retry transaction. Every physical molecular macrostep and complete half interval
receives the original endpoint entropy gate. Its complete pre-step source
account retains physical chemical entropy constants; only the reference effect
of bounded numerical source-continuity defects is removed.

One actual RHS budget spans both halves and discarded compositions. Exhaustion
is terminal. All fields, clocks, finite-port receipts and operator reports
publish only after complete temperature, balance, entropy and final cancellation
checks. `StabilizedReactiveNetworkEntropyReport` retains the complete entropy
result, both accepted RKL2 intervals and total actual work. Native source and
transport owners retain their existing workspaces and original caloric laws.

This entry admits stationary closed constant-area geometry with common species
support at active molecular faces. Baths, moving volumes, other sources and
adaptive entropy composition remain separate paths. The physical macrostep cap
is an explicit accuracy choice. Linear RKL stability and symmetric splitting do
not establish nonlinear positivity, a fully discrete coupled entropy theorem
or global temporal order. Focused G1/G4 regressions check complete variable-Cp
chemical entropy, heat-release retry budgets, final cancellation and reuse.

### Complete source entropy in open reacting networks

`ReactingNetwork1d::step_with_reservoirs_and_entropy` admits the original
C(dt/2)-A(dt)-C(dt/2) transaction with stationary prescribed pressure reservoirs.
The native characteristic closure, chemistry, engineering flux and bounded
complete-duration retries are retained. Every branch cell and finite plenum
contributes actual source entropy. All-None reservoirs retain the closed path.

The complete account is anchored before chemistry. Each accepted physical
exterior pair retains its actual SSP species/momentum/referenced-energy
transfers, signed time-area measure and absolute stage species throughput.
These inventories are contracted against each original reservoir effort in
the immutable pre-chemistry basis; rounded component entropy reports are not
summed. Internal finite-plenum ports and periodic faces supply no exterior
entropy. Advected enthalpy already contains pressure-flow work. Outflow uses
the receiving reservoir's thermodynamics and includes coupling irreversibility.

Actual chemical entropy constants remain physical. Only bounded numerical
species-source continuity defects lose their reference contribution. Zero
reservoir activity refuses nonzero transfer at either SSP stage. The complete
endpoint gate and final cancellation poll precede publishing fields, clocks
and receipts. G1 tests reconstruct variable-Cp state entropy and reservoir
efforts independently, retain a nonzero chemical entropy datum, admit entropy-
decreasing outflow with positive generation, and exercise hot-reaction retries.
G4 checks final cancellation, zero-activity refusal and reuse.

This entry requires fixed geometry and excludes molecular baths, moving
volumes and additional wall/valve sources. It measures finite-step entropy
admission of the existing engineering model; no fully discrete stability or
global temporal-order theorem is claimed.

## Molecular N2/O2 shockwave manifold consumer

`cargo run -p fs-gas --features experimental-euler --example mixture_manifold -- 64 molecular`
evolves the pressure/composition pulse with `DissipatingNetwork1d<KineticTransport>`.
The existing default remains inviscid and retains its conical daughter. Molecular
mode explicitly selects three stationary cylinders of lengths 0.6/0.5/0.7 m and
radii 0.010/0.008/0.008 m, a 30 mL rigid plenum, and a declared 0.010 m plenum-side
mixing length at each port; each gradient includes the endpoint half-cell width.

Calorics, masses and reference energies remain the original airNASA9 source bank.
The dilute nonpolar kinetic law uses original Cantera v3.2.0 `data/gri30.yaml`
collision data (blob `9a62004ec17bbec8ac11e97785932d426dd58610`): N2 sigma
3.621 angstrom and epsilon/k 97.53 K; O2 sigma 3.458 angstrom and epsilon/k 107.4 K.
It declares Wilke viscosity, estimated Eucken conductivity, zero bulk viscosity,
200--1000 K and 1000--1e6 Pa. This does not resolve radial wall friction/heat
transfer, dense gas, Soret/Dufour or rotational relaxation, and the compact
mixing distance is a constitutive input rather than a resolved plenum flow.

Fresh native runs with Maxwell--Stefan evaluator version 2 completed 4 ms:

| Cells per branch | Wave steps | Molecular substeps | Maximum species residual (kg) | Physical-energy residual (J) |
| --- | --- | --- | --- | --- |
| 16 | 168 | 336 | 1.36e-20 | +5.19e-15 |
| 64 | 669 | 1338 | 1.09e-19 | +1.43e-14 |

Both runs had zero rejected molecular substeps, coupling retries, HLL fallbacks
and timestep halvings. Every sampled field was finite, with positive pressure,
source-admitted temperature, and species fractions in [0,1] summing to one at
printed precision. The sampled minimum N2 fraction was 0.765013831926 at 16 cells
and 0.766338161781 at 64 cells, below initial air's 0.76708: these are not
composition maximum-principle, entropy, convergence-rate or experimental
accuracy certificates. Both daughters receive pressure/composition waves and
nonzero conservative molecular port transfers. The rebuilt default 16-cell
inviscid run retained its prior 149 steps and -3.17e-14 J physical-energy residual.

## Stationary wall drag and finite wall heat storage

`multicomponent::wall_exchange` applies a source interval to one original
mixture cell without changing species or the flow clock. Linear drag uses its
exact frozen-rate exponential; Darcy drag uses its exact frozen-factor
quadratic decay with an explicit hydraulic diameter. The gas retains total
energy while lost kinetic energy becomes sensible energy, and the stationary
support receives the opposite signed axial impulse. No second viscous-heating
term is added to total energy.

Opt-in `WallDrag::MovingLinear` and `WallDrag::MovingDarcy` prescribe constant
axial wall velocity U over a source interval at fixed gas volume. They exactly
relax relative velocity u-U using the same frozen linear/Darcy laws. Gas receives
mechanical work U times its actual accepted momentum change; the wall/drive
receives the opposite impulse and work. Loss of relative-slip kinetic energy
becomes sensible gas energy once, even during entrainment or braking.
`mechanical_work_to_gas_j` remains separate from `heat_to_gas_j` and
`drag_heating_j`, and frozen branch reports retain `wall_work` separately from
`wall_heat`. Both frozen and reacting/molecular mixture network balances include
drive work. Original caloric wall storage and material-derived coefficients use
the same owner; coefficient refresh preserves wall motion and admits Reynolds
and Mach numbers in the wall frame. Stationary cases retain their original
arithmetic. This is a prescribed sliding-wall source, not resolved wall inertia,
moving geometry or an unsteady boundary-layer model; zero-momentum plenums
continue to refuse drag. Four numerical-kernel and eight native integration
regressions cover analytic slip relaxation, frame covariance, weak/strong drag,
finite calorics, network composition and transactional refusal.

`WallDrag::LinearQuadratic` applies authored linear and inertial resistance
simultaneously to slip s = u-U: ds/dt = -a*s - b*abs(s)*s, with nonnegative
a [1/s], b [1/m] and a prescribed constant wall velocity U. Its single exact
frozen-coefficient update is s(t) = s(0)*exp(-a*t) /
[1 + b*abs(s(0))*(1-exp(-a*t))/a], using the continuous a=0 limit. Stable
exponential differences retain weak impulses; scaled rational decay admits
strong resistance without splitting the two forces or clipping the residual
slip. The existing accepted-impulse work and slip-heating account is shared.
The zero-coefficient limits are the original linear and quadratic laws.
G1 regressions compare an independent reciprocal-slip solution, subdivision,
frame transformations and weak/strong limits; the mixture-network source test
also checks actual species, temperature, support impulse and drive-work balances.
No coefficient, transition law or applicability range is inferred.

`exchange_cell_with_recoil` couples one original mixture cell to a finite
`SlidingSupport` on the same axis at fixed gas volume and species mass. The
existing exact linear/quadratic slip map uses the reduced mass m*M/(m+M), so
both bodies respond to the same interval. Accepted gas impulse I gives the
support -I, wall work I*(U0-I/(2*M)), and nonnegative slip loss
-I*(s0+I/(2*reduced_mass)). An authored `gas_heat_fraction` partitions this
loss; a nonzero wall share requires finite `WallThermal::Lumped` storage.
The original caloric heat owner then advances the contact without resetting a
flow clock. Returned gas/support endpoints and actual finite-wall storage must
pass independent momentum and complete energy budgets before return.
`RecoilUpdate` exposes both slip-heat shares; the original `WallCellReport`
retains total slip loss and accounts for the wall share as heat leaving gas,
so its sensible-energy identity remains additive. Inputs are immutable on
refusal or cancellation. This owner binds one support to one gas amount;
shared multi-cell inertia and moving geometry require their own coupled owner.
G1 regressions cover reduced-mass analytic decay, subdivision, large mass ratios,
frame covariance, heat fractions 0/0.5/1 with subsequent conduction, actual
combined storage, and late cancellation. Friction coefficients and the heat
fraction are authored; the thermal step retains its original first-order law.

`exchange_cell_with_recoil_and_heat_storage` retains the same mechanical owner
and accepts the support's original `MaterialWallCaloric`. Frictional wall heat
is deposited through that source's exact capacity integral or native energy
coordinate before the existing implicit heat exchange. Temperature-dependent
capacity and latent internal-energy plateaus therefore remain active during
support recoil. The energy audit uses the complete original-to-returned wall
storage difference, not a frozen capacity or a sum of rounded stage changes.
The source must explicitly declare internal energy or negligible solid
expansion; unqualified enthalpy and sensible-capacity storage cannot close this
mechanical balance. Generic heat-only wall exchange still accepts those source
coordinates and retains their declared interpretation. The thermal amount is
already part of the authored support inertia and is not added to it again.
Original source-domain exits and cancellation return neither endpoint. Focused
regressions cover analytic variable-capacity slip heating, variable-Cv mixture
exchange, a latent plateau, the constant-capacity limit, and atomic refusals.

`SharedRecoilHeatWorkspace` couples several original gas amounts to one finite
support and one original `MaterialWallCaloric` inventory. Each contact supplies
its frozen resistance, gas/wall frictional heat fraction and conductance. The
existing shared mechanical solver advances one common inertia; the combined
wall share is deposited once through original nonlinear or latent storage.
The existing shared-wall scalar solve then exchanges heat simultaneously with
every original gas composition against that one wall endpoint. Whole returned
gas, support kinetic and material energy changes must close the original
budget. The result retains separate mechanical and thermal-stage reports,
per-cell work, heat and wall shares, and the admitted storage interpretation.
Allocation occurs at preparation; bounded stepping borrows all gas species and
publishes no caller state. Cancellation or refusal permits exact retry from
the supplied gas, support and wall inputs. The same internal-energy or explicit
negligible-expansion requirement applies. Frozen coefficients, fixed geometry
and first-order thermal splitting remain explicit; network bindings continue
to own their separate source transactions.

Thermal contacts may be adiabatic, a prescribed-temperature reservoir, or a
finite lumped wall with explicit heat capacity in J/K. A contact's conductance
is the complete exposed h*A in W/K. A safeguarded backward-Euler temperature
solve evaluates the original mixture's NASA-9 energy and Cv rather than a
frozen-gamma substitute. Finite-wall storage is eliminated in the scalar solve
and its actual represented temperature change is included in the energy audit.
Reservoir heat is external; gas/finite-wall exchange is internal. The step
refuses source-domain or residual-budget violations rather than clipping gas
temperature or creating heat.

The common wall caloric path locally refines gas temperature when the flow
inverse is too coarse for the heat equation. It uses the original energy and
Cv queries with the actual or predicted heat scale and `G*dt/Cv` sensitivity,
both before forming a passive bracket and at the returned thermal endpoint.
Original source domains, energy references, conservative states and caller EOS
controls are unchanged. The wall's own iteration limit bounds refinement;
shared exchanges charge these queries to their existing nested source budget.
Zero conductance and exact caloric states need no additional inversion.
Native h/u heat bounds account for the represented mass conversion before a
trial source query: only the proposed extensive bracket contracts toward zero
until its reconstructed state lies in the declared source interval. Accepted
stored energy is neither clipped nor reprojected.

`network::wall_exchange::NetworkWallOperator` binds one distinct contact to
every branch cell and finite plenum. It stages all momentum/energy and wall
temperature updates in prepared buffers, preserves each species exactly, and
publishes the entire source interval only after its final cancellation check.
Failure leaves accepted gas, walls, clocks and reports intact. Its accepted
support impulses remain per branch axis. A successful standalone source update
clears an old advection report that no longer describes the changed fields.
The source operator accepts the existing stationary mesh geometry; the separate
molecular operator retains its constant-area restriction.

The network has a shared finite-wall trial workspace for complete coupled
retries. Candidate walls start from accepted storage once per trial, retain
both wall half steps, and publish only when the composing owner accepts every
stage. The report separates heat into gas, actual finite-wall energy change,
external reservoir heat, kinetic energy thermalized, and numerical residuals.

These are Estimated lumped/averaged closures with explicitly supplied
coefficients. They do not resolve developing, turbulent or unsteady shock
boundary layers; solid conduction between contacts; shared multi-cell walls;
radiation; moving-wall work; or a thermodynamic/entropy theorem. Backward-Euler
heat integration is first order even when surrounded by symmetric source
composition. No throughput or experimental-accuracy claim is made.

Focused Rust regressions accompany the source/network implementation for an
independent two-body heat solution and timestep refinement, exact linear/Darcy
drag and support impulses, reservoir exchange/source refusal, and late
cancellation followed by reuse. Native compilation and execution are blocked in
this session: the execution environment has no usable connection. Consequently
these new regressions have not been run, and no DSR, Cargo, native-suite or
physical-consumer pass is claimed for this extension.

### Source-coupled viscous slip relaxation

`exchange_cell_with_source_drag` advances the original mixture's mechanical
wall source with coefficients reevaluated at its evolving conservative state.
The borrowed callback supplies an admitted `WallDrag`; it does not supply a
replacement caloric model. Gas amount, composition, volume and prescribed wall
velocity remain fixed during this interval. At every candidate state the
original gas owner admits its density, temperature and pressure. No wall heat
transfer is performed by this mechanical kernel.

For gas momentum density p, wall speed U and slip s = p/rho - U, the invariant
is E - U*p. Each accepted impulse therefore determines drive work and physical
slip dissipation from the actual returned endpoints. Stationary walls retain
total gas energy exactly; moving walls transfer U times the gas impulse.
Relative kinetic energy lost becomes gas sensible energy once. There is no
extra heating source and no conversion of a numerical error estimate to heat.

A source evaluation at an exact frozen-law half step predicts the midpoint
coefficients; the original linear/Darcy/combined exact map then advances the
entering state for the complete candidate interval. One full candidate is
compared with two half candidates, and the refined state is accepted only
within independent dimensional slip and physical-heating budgets. The
full-versus-refined difference divided by three estimates the smooth
second-order local endpoint defect. Work is bounded and cancellation/source
refusals publish no endpoint.

A separate stage-order check requires the predicted midpoint to lie between
the entering and returned slip magnitudes. This rejects the false agreement
that can occur when both coarse and refined frozen predictors cross a sharp
viscosity decrease. Constant positive coefficients retain their exact stiff
decay. The embedded defect and stage-order check are numerical acceptance
tests, not a certified error bound for arbitrary unresolved constitutive
features or a global thermodynamic theorem.

The returned source receipt separates accepted/rejected work and error
estimates from the original physical wall report. Focused G1/G3 regressions
use an independent closed-form temperature-proportional viscosity solution,
constant-coefficient limits, opposite slip and moving-frame transformations,
smooth sharp-decrease refusal/refinement, source-domain exit, bounded work
and cancellation. Thermal-contact splitting and radial boundary-layer memory
remain separate owners.

### Live material resistance with independent support recoil

`exchange_cell_with_source_recoil` shares the source-coupled scalar driver
with prescribed-wall drag and retains the original independent recoil owner.
For gas mass m, support mass M, reduced mass mr = m*M/(m+M), entering slip s0
and actual gas impulse I, the original mechanics owner determines:

- support velocity U = U0 - I/M;
- gas work W = I*(U0 - I/(2*M));
- total physical slip loss D = -I*(s0 + I/(2*mr));
- gas energy E = E0 + (W - (D - f*D))/V for the fixed authored gas fraction f.

Every trial is anchored to the original gas/support state through that same
accepted-impulse calculation. The original gas primitive and candidate support
are passed to the constitutive callback; its returned drag must use that
support velocity. The callback cannot change support mass or heat fraction.
Slip defects include both bodies' response, and heating estimates control
total reduced-mass loss before partition, including an all-wall heat share.
Actual represented support momentum and complete energy remain independently
audited; rounding residuals are never repaired by an impulse or heat source.

The original recoil completion deposits D - f*D through the authored constant
capacity or original `MaterialWallCaloric` integral/native energy coordinate,
then runs the original thermal-contact owner. Internal-energy or explicitly
negligible-expansion admission remains required for material support storage.
The thermal material amount belongs to the authored inertia and is not added
to it again. The final gas/support state is admitted by the original resolver
after thermal completion, with that query charged to the same whole-call
source budget. No endpoint is returned after refusal or cancellation.

The network binding admits independent supports in either binding order and
uses their existing trial/accepted `SupportExchange` accounting. Transport is
queried in the actual moving-support frame; after mechanical work and friction
heat partition, the common source seam refreshes the local gas film before
ordinary caloric/conduction exchange. The original admitted gas primitive is
reused in each constitutive stage, avoiding a duplicate caloric inversion.
Shared finite inertia still requires its collective vector owner and is
refused on selected source-coupled branches. Shared thermal inventories retain
their original separate conductances and heat solve.

G1/G3/G4 regressions compare an independent temperature-proportional-viscosity
reduced-mass solution over light/moderate/heavy supports, opposing slip, frame
boosts and gas heat fractions 0/0.5/1. Constant linear/Darcy/combined resistance
recovers the original recoil result. Original nonlinear wall storage,
post-partition gas-film refresh, final thermal source refusal, late cancellation
and deterministic complete-network retry are covered. The source controller
does not adapt the subsequent thermal-contact equation or claim a global
second-order coupled thermal method.

## Wall sources composed with molecular and reacting networks

### Current material state in developed pipe-wall losses

`NetworkWallOperator::with_circular_laminar_transport` binds the existing
`TransportLaw` (including source-card `MaterialTransport` and original kinetic
transport) to selected circular branches. Full caloric/species-source identity
must match the network. At every wall source half it obtains temperature,
pressure, density and composition from the actual conservative cell, reuses the
original molecular mole-fraction closure, and queries the same source-owned
viscosity and conductivity. Thus new reaction products, changed inlet mixtures,
pressure or thermal evolution change wall drag and cooling without a caller
re-entering property scalars or a separate reacting-wall solver.

The declared `CircularLaminarWall` supplies diameter, Nusselt number, supported
Reynolds interval (within 0..2300), and maximum Mach (at most 0.3). Diameter and
cell volumes must agree with the original constant-area circular mesh. Drag is
`32 mu/(rho D²)`; complete cell conductance is `Nu k (4 V/D)/D`, where `4 V/D`
is that cell's wetted area. The Nusselt number and fully developed applicability
are explicit engineering inputs, not values inferred from a gas name. This
closure does not establish entrance length, negligible axial conduction,
high-frequency oscillatory boundary layers, turbulent flow or a measured pipe
correlation. Source coefficients are evaluated at bulk gas temperature and, by default,
frozen during each source half; that constitutive lag remains first order.
The opt-in source-coupled drag binding below resolves the mechanical feedback
within each interval while retaining the existing separate thermal owner.

Both source-step endpoints must satisfy the original property law and declared
Reynolds/Mach limits. A changed product outside its material-card domain rejects
the complete trial, including earlier cells, wall temperatures, reaction updates
and reports. The existing exponential drag, implicit finite-wall thermal solve,
species conservation, support impulse and energy accounting remain unchanged.
No molecular or chemical heat source is added twice. Finite wall heat capacity
is still an explicitly authored constant; this binding updates gas-side
transport, not solid-wall calorics. Unbound branches and stationary plenums
retain the existing explicitly authored contacts. Original law workspaces may
change during refused trials, but carry no accepted physical state.

Focused G1/G3 regressions check independent drag/finite-wall heat equations,
density and transport-data substitution, exact caloric-source and geometry
admission, heated-endpoint rollback, and positive reaction-product refresh
through `ReactingNetwork1d::step_with_wall`. An original material-card refusal
after reaction also checks whole-network/finite-wall rollback. These are
manufactured numerical integration tests, not experimental combustion or
convection validation.

### Viscosity feedback within developed wall-source intervals

`NetworkWallOperator::with_source_coupled_drag` opts selected developed-flow
branches into the original source-coupled mechanical kernel. Bind it after
the original `TransportLaw` and geometry, before evolution. At every admitted
candidate momentum/energy state, the binding evaluates that same law at the
actual original-mixture temperature, pressure and composition, then recomputes
the existing geometry- and Reynolds-dependent resistance. Circular laminar,
rough/turbulent, Gnielinski, Churchill, rectangular and annular models retain
their original correlations and declared domains; no viscosity fit, source
table clone, inferred transition or clipping is introduced.

After friction changes the gas state, the binding refreshes the local film
conductance before handing off to the original finite/reservoir/caloric wall
and solid-conduction owners. Shared thermal inventories retain their authored
port conductances and original collective heat solve. Species, branch
impulses, physical slip heat, drive work, finite-wall energy, cancellation and
whole-network retry remain on the existing transaction. Adaptive source work
limits apply independently to each selected cell and source interval.

The binding admits stationary walls, prescribed motion and independent finite
supports with their original frictional heat partition and material storage.
A shared inertial support on a selected branch is refused in either binding
order; unselected branches retain their existing support owners. The scalar
mechanical solve is second order for smooth live coefficients, while the
separate backward-Euler thermal contact and outer coupled system retain their
existing order limitations. Conductivity is refreshed after mechanical
heating but is not solved jointly with temperature during the thermal step.
These remain fully developed averaged wall closures, without unsteady radial
profile storage or a resolved turbulent boundary layer.

### Rectangular developed-flow wall resistance

`NetworkWallOperator::with_rectangular_laminar_transport` binds a prepared
`RectangularLaminarWall` to selected constant-area rectangular branches.
The original widths, heights, Nusselt inputs and applicability are available
through `rectangular_walls`. Circular entry points and their arithmetic remain
unchanged. Both shapes share the source-owned transport query, original wall
transaction, exponential relative-slip drag, moving-wall/finite-support work,
finite thermal storage and atomic publication.
`with_developed_transport` combines circular and rectangular branch models
under one original transport law. `developed_walls` retains the complete mixed
axis; the circular/rectangular slice inspectors apply to homogeneous bindings.

For smaller side h, larger side w and a=h/w, the no-slip Fourier solution gives
`Q = w h^3 C(a) dp/(12 mu L)` and hence bulk drag rate `12 mu/(rho h² C(a))`.
The complete side-wall correction is retained; it is not circular Poiseuille
with a hydraulic-diameter substitution. A fixed six-term exponential-difference
sum evaluates the shape factor once at preparation; the omitted correction is
below 1e-23 on 0<a<=1. Source intervals refresh current gas viscosity and density
without repeating that series. Reference: Mortensen, Okkels and Bruus (2005),
Phys. Rev. E 71, 057301, equations (3), (10), (11),
<https://doi.org/10.1103/PhysRevE.71.057301>.

Actual area w*h and integrated cell volumes must agree with the original mesh.
Hydraulic diameter `Dh=2wh/(w+h)` controls Reynolds admission, authored Nu and
cell wetted area `4V/Dh`. Nu must apply to the declared rectangular aspect ratio
and thermal wall conditions; no circular heat-transfer relation is inferred.
All four walls share their contact temperature and axial velocity. Reynolds
support lies within 0..2300 and wall-relative Mach is at most 0.3. These bounds
do not establish laminarity, entrance development or quasi-steady validity.
Slip, rarefaction, turbulent rectangular flow and unsteady transverse profiles
remain outside this model. Source coefficients retain the existing first-order
interval lag and endpoint domain checks.

G1/G3 tests compare square-duct resistance, an independent long Fourier sum,
parallel-plate convergence, swapped sides, signed flow and moving-wall boosts.
They also check actual perimeter heat exchange, density/material substitution
and kinetic-to-thermal conversion without a second heat source. G0/G4 checks
cover source/geometry/regime refusal and cancellation with complete rollback.

### Concentric annular gas-wall resistance

`AnnularLaminarWall` and `with_annular_laminar_transport` extend the original
gas-wall owner to concentric passages, including narrow cooling/seal gaps.
`DevelopedWallModel::Annular` shares a transport source and one conservative
transaction with circular and rectangular branches. Actual area is
`pi (R-r)(R+r)`, hydraulic diameter is `2(R-r)`, and both wetted surfaces enter
the cell conductance. Original source viscosity and density refresh drag;
the explicitly authored Nu retains the annular radius-ratio and thermal-wall
assumptions. Both cylinders share one contact temperature and axial velocity.

The no-slip axial solution gives
`Q = pi dp [R^4-r^4-(R²-r²)²/ln(R/r)]/(8 mu L)`.
Preparing its exact mean-flow resistance once avoids a logarithm or series at
source-step time. For r/R>=1/2, a fixed positive series in `(R-r)/(R+r)` avoids
catastrophic cancellation and preserves `12 mu/[rho (R-r)²]` in the thin-gap
limit. A nonzero vanishing core retains logarithmic drag even when r/R
underflows; only exactly zero inner radius selects circular arithmetic.
Reference: Huilgol and Georgiou (2020), Newtonian specialization of equations
(18), (22), (23), <https://doi.org/10.1016/j.jnnfm.2020.104401>.

The original interval-frozen source model, Reynolds ceiling 2300, wall-relative
Mach ceiling 0.3, endpoint domain checks, heat/work ownership and rollback
remain unchanged. Eccentricity, differential cylinder velocities, separate
inner/outer wall temperatures, rotation, slip, turbulent annular flow and
resolved unsteady transverse profiles are outside this closure.

G1/G3 kernel tests check independent 100-digit integrated-velocity values,
thin-gap resistance, circular-limit bit agreement, logarithmic small-core drag
and material/length scaling. All six kernel tests were executed by compiling
the actual source with the production `fs-math` library. Network regressions
for both-wall heat transfer, mixed shapes, boosted wall motion, energy accounting
and late regime refusal are authored but await native integration execution;
the isolated kernel result is not evidence that those regressions passed.

### Current material state and roughness in turbulent pipe-wall friction

`NetworkWallOperator::with_circular_transport` accepts one explicit
`CircularWallModel` per selected branch. Laminar and turbulent branches use the
same original `TransportLaw`, conservative gas fields, thermal storage,
finite-support mechanics and atomic source transaction. The laminar convenience
constructor retains its existing Poiseuille law; `circular_walls` exposes the
complete selected models in branch order.

`CircularTurbulentWall` resolves the Darcy factor from the implicit
Colebrook-White equation using actual density, source-owned viscosity and
wall-relative velocity. It requires authored equivalent sand-grain roughness,
the original circular diameter, an explicit Reynolds interval within
4000..1e8, roughness/diameter in 0..0.05 and a Mach ceiling at most 0.3.
Equivalent sand roughness is not inferred from a bulk material name or equated
to a surface's RMS height. Its Nusselt number remains an independent declared
heat-transfer input. Reference: Colebrook (1939),
[DOI 10.1680/ijoti.1939.13150](https://doi.org/10.1680/ijoti.1939.13150), and
[NASA KSC-STD-Z-0017 Rev A](https://standards.nasa.gov/sites/default/files/standards/KSC/A/0/KSC-STD-Z-0017_RevA_Final.pdf),
section 7.4.1, equation 4.

The allocation-free scalar solve uses a strictly monotone reciprocal-square-root
residual, deterministic logarithms, safeguarded Newton updates and a 32-iteration
cap. The existing exact frozen Darcy step then advances momentum. Prescribed
wall work, opposite support impulse, slip heating and wall heat retain their
original separate energy accounts. The factor refreshes at every source
interval; entering and leaving material/flow states must both be admissible.
A domain exit, source refusal or cancellation publishes no partial gas or wall
state. The Colebrook option supplies no interpolation through transitional flow.

This is an Estimated developed-flow correlation with first-order constitutive
lag, not resolved turbulence, an unsteady boundary-layer model or a
high-compressibility correction. G1/G3 regressions check the original Colebrook
equation, smooth/rough limits, viscosity/roughness changes, reverse slip,
moving-frame work and heat, mixed laminar/turbulent branches and refreshed
factors. An 8/16/32-step ladder compares native braking with an independent RK4
solution of the continuously refreshed Colebrook law and checks first-order
constitutive convergence. A G4 endpoint-domain regression checks the actual cells, finite walls
and previous report after refusal.

### Turbulent gas-wall heat transfer from the current mixture

`CircularWallModel::Gnielinski(CircularGnielinskiWall)` explicitly selects
coupled smooth-pipe friction and turbulent heat transfer. Each original source
interval evaluates `Re = rho |u-u_wall| D / mu` and `Pr = mu Cp / k` from its
current mixture. `Cp` comes from the original NASA calorics and composition;
it is neither `Cv` nor a constant-gamma substitute. The smooth-pipe Darcy
factor `(0.79 ln(Re)-1.64)^-2` drives both the momentum loss and Gnielinski's
Nusselt relation. The existing gas/finite-wall transaction then uses the
cell's complete wetted area `4 V/D`, including its original solid resistance
when material wall conduction is bound. No new heat source duplicates the
existing conversion of relative kinetic energy to sensible energy.

Reference: Gnielinski (1976), *International Chemical Engineering* 16(2),
359–368, turbulent pipe-flow heat-transfer correlation.
The implemented range is deliberately narrower: authored Re within 1e4..1e6,
Pr within 0.5..200, and wall-relative Mach at most 0.3. Original circular
geometry must agree and branch length/diameter must be at least ten. Fully
developed smooth flow and modest wall-to-bulk property variation are modeling
assumptions; the length check does not prove development. This path does not
extend the heat correlation to rough pipes, transition, swirl, strong property
variation or oscillatory boundary layers. The existing Colebrook option keeps
its independently authored Nusselt number.

G1 regressions check the published friction/heat relations and the Pr=1 limit.
G3/G4 native-network tests check live-property and slip response, finite-wall
energy exchange, moving-frame consistency and complete rollback on domain
exit. The constitutive coefficients retain the source owner's first-order
interval lag. These are numerical checks of an Estimated correlation, not
experimental gas-pipe validation or a resolved-turbulence claim.

### Explicit pipe friction through startup, transition and stopping

`CircularWallModel::Churchill(CircularChurchillWall)` opts a branch into
Churchill's continuous engineering friction relation. Its current original
mixture density, viscosity and wall-relative speed set Reynolds number; the
caller supplies effective sand roughness and an independent Nusselt number.
Authored Re must lie within 0..1e8, roughness/diameter within 0..0.05, and the
Mach ceiling must not exceed 0.3. Original circular geometry is checked.
The restricted laminar, Colebrook and Gnielinski models retain their domains.

Reference: Churchill (1977), *Friction-factor equation spans all fluid-flow
regimes*, equations 5 and 18. The paper's factor is `f_D/8`. The allocation-free
deterministic kernel evaluates the finite product `f_D Re`, with rest limit 64,
and retains small excesses above Poiseuille resistance without cancellation.
The existing exact linear-quadratic source keeps `32 mu/(rho D²)` as its linear
term and freezes only its coefficients over each interval. It uses the original
conservative wall mechanics, heat, finite-support and rollback owners.

This is an Estimated quasi-steady developed-pipe approximation with first-order
constitutive lag. Transition values are uncertain and do not resolve disturbance
dependence, intermittency, hysteresis or unsteady boundary-layer memory. Friction
does not supply a heat-transfer blend: the authored Nusselt number and its
applicability remain independent. Focused regressions exercise the published
equation, rest and rough limits, native braking across all three regimes with
time refinement, live properties, moving-wall energy, and late source refusal.
These checks do not establish experimental pipe-flow validation.

### Existing source composition

`DissipatingNetwork1d::step_with_wall` performs W/2, D/2, A, D/2, W/2.
`ReactingNetwork1d::step_with_wall` performs W/2, C/2, A, C/2, W/2;
`step_with_dissipation_and_wall` surrounds the existing D/2, C/2, A, C/2, D/2
composition with both wall halves. The existing methods continue without W.
Every half covers the actual accepted wave duration. A shorter advection step
or a refused source stage resets gas fields and finite wall storage before
retrying the complete composition.

Both wall heat transfers enter the gas boundary-energy balance. The wall report
separately retains actual finite-wall storage changes and external reservoir
heat, so a gas-plus-wall audit can distinguish internal exchange from external
energy supply. Stationary-wall kinetic thermalization is diagnostic only and
does not enter the total-energy source a second time. Final branch and plenum
temperatures must still fit the chemistry mechanism's admitted interval after
the final wall half. Gas, molecular reports, wall temperatures and both wall
halves' aggregate support impulses publish only after the final cancellation
check.

The manifold example accepts `molecular-wall`:

```bash
cargo run -p fs-gas --features experimental-euler --example mixture_manifold -- 64 molecular-wall
```

It retains the original N2/O2 caloric and kinetic sources, and adds explicitly
illustrative fixed inputs: Darcy factor 0.02, h=50 W/(m² K), wall areal capacity
2000 J/(m² K), initial wall temperature 300 K, cylindrical sidewall areas and an
authored 0.005 m² plenum wall area. It independently integrates final retained
wall temperatures to report a complete gas-plus-wall energy balance alongside
gas-only and molecular-port budgets. These coefficients are not a named
material calibration or an unsteady shock-layer model.

Focused coupled regressions exercise closed gas-plus-wall storage, zero-wall
identity, exact half-duration accounting and late cancellation/reuse. Additional
cases check an independently integrated nonconstant mixture Cv polynomial and
simultaneous reaction, molecular transport and finite-wall energy storage.
These new Rust tests and the new example mode are committed but have not been compiled
or executed in this session: `dsr quality --tool frankensim` could not launch
because the execution environment is unavailable. No native/DSR pass, measured
wall attenuation, global second-order accuracy or full entropy certificate is
claimed.

## Finite-rate chemistry in common-pressure gas layers

`multicomponent::compartment::kinetics::ZoneKinetics` binds the complete original
gas/kinetic source and a prepared room geometry. Each selected present layer
uses the existing `IsochoricReactor` at its actual volume and physical internal
energy, including formation terms. The resulting time-integrated species
changes become canonical `LayerRates`; their only energy term converts the
original numerical species references. No second combustion-heat source is
added. Both layers then use the original simultaneous pressure/interface-work
solve. The enclosing transport owner advances the physical clock once.

The prepared reactor and room candidate are reused. Actual rate reconstruction,
mass and physical-energy residuals are checked before publication. Gas/source
and geometry mismatches, exhausted work, domain failures and late cancellation
leave the accepted room unchanged. The kinetic temperature interval is checked
again after mechanical relaxation: compression can invalidate a layer even
when its homogeneous reactor endpoint was admissible. A shorter retry starts
from the supplied accepted state, independent of previous workspace contents.

The chemistry report retains original per-layer entropy, conservation, local
error and adaptive-work diagnostics. The selected fixed-volume chemistry
method is followed by first-order mechanical/source composition; the homogeneous entropy check is not a
fully discrete theorem for the moving-interface or ventilated system. Supplied
mechanisms retain their own source and validity limitations. No ignition,
flame-structure, CO/soot-yield or experimental accuracy claim is inferred.

Eight new regressions cover analytic finite-time reversible relaxation,
unequal reacting layers with variable heat capacities and actual interface
work, disabled/absent layers, trial replay, formation-reference invariance,
source/domain/work refusals, post-relaxation domain rejection and late
cancellation. All 23 focused native compartment tests passed; that is separate
from the unavailable full-workspace DSR gate.

## Reduced finite-rate Newton system and independent hydrogen reactors

The original extent integrator selects the smaller of the full reaction
matrix and an equivalent species-plus-temperature Schur system. With ten
species and 29 reactions it factors 11 rather than 29 unknowns. The additional
thermal coordinate retains the original compensated temperature derivative
with respect to each reaction extent, preserving elemental energy-reference
invariance. Every reaction extent is reconstructed, and bounded residual
correction checks the solution against the original full linear equations.
The original pivoting, positive-amount, caloric-domain, nonlinear, energy,
element/charge, entropy, local-error, retry and cancellation checks still apply.
No measured wall-clock speedup or higher integration order is implied.

Five focused reduced-system cases compare independent full/reduced equations,
stiff cancellation recovery, exothermic conservation, elemental datum shifts,
and deterministic workspace reuse. The original gas library passes 151 tests.
A further independent Cantera 3.2.0 test uses the complete sourced hydrogen
mechanism, continuous calorics, no radical seed, 1 litre at 101325 Pa and
H2:O2:N2 = 1:1:3.76. Initial 1000 K advanced 1 ms and 1200 K advanced 0.1 ms
agree with Cantera within 3 K, 0.2% pressure and 0.5% + 1e-12 mol per species;
physical internal energy agrees within 2e-8 J and atom/mass inventories close.
For the 1200 K case, tightening local relative tolerance from 2e-4 to 2e-6
reduces the combined temperature/species endpoint error by over 70%.
`data/reference/thermochem/hydrogen_cantera32.py` reproduces the independent
reference using a separately installed, pinned development-only Cantera.
These are numerical source/integration comparisons, not experimental ignition
validation, a fully discrete coupled entropy theorem, or a second-order claim.

## Optional second-order L-stable homogeneous chemistry

`IsochoricKineticsSolve::method` explicitly selects
`IsochoricKineticsMethod::BackwardEuler` (unchanged default) or `Sdirk2`.
The latter implements Alexander's two-stage stiffly accurate method with
gamma = 1 - 1/sqrt(2), A = [[gamma, 0], [1-gamma, gamma]], b = A[1].
Reference: R. Alexander, SIAM J. Numer. Anal. 14(6), 1006–1021 (1977),
https://doi.org/10.1137/0714068. This is a second-order L-stable method for the
smooth homogeneous reactor ODE, not an order claim for any surrounding split
transport, radiation, pyrolysis or pressure-interface composition.

Both stages reuse the original constrained extent Newton owner, source fixed
U,V, smaller factorization, caloric inversion, positive-amount and work gates.
The second stage keeps its known reaction forcing separate from the physical
baseline and starts Newton from the admitted first-stage state; an explicit
shifted baseline may be nonphysical and is never used as a state. Each stage
must pass an entropy comparison with its physical step baseline. Stage two is
not incorrectly required to exceed stage one's entropy. These are actual
acceptance gates, not unconditional Runge–Kutta positivity/entropy theorems.

One reaction-sized forcing vector is reserved only for SDIRK2 at construction;
advances allocate no new stage storage. The two-half-step result is retained
without Richardson extrapolation, with the error difference divided by three
and a deterministic cube-root step controller. Backward-Euler behavior and its
square-root controller remain the default.

Four new native tests demonstrate fixed-step second-order convergence against
analytic reversible relaxation and its independent stability function, stiff
damping, physical energy/elemental datum invariance, cancellation with replay,
and two independent Cantera 3.2.0 hydrogen endpoints. All 156 native gas tests
pass. At identical local relative tolerance 1e-5, the 1200 K/0.1 ms hydrogen
fixture uses 322 accepted SDIRK2 steps and 3714 Newton iterations versus
4485 and 19983 for backward Euler; temperature errors against Cantera are
0.01556 K and 0.04168 K respectively. These are fixture-specific work counts,
not a universal speed or accuracy guarantee.

## Constant first-order half-life evolution

`IsochoricKineticsMethod::ExactFirstOrder { limits }` explicitly selects the
shared positive `fs_la::decay::reaction::ReactionPlan` for a constant chemical
generator. Every source channel must be irreversible, have exactly one
order-one reactant, and use an Arrhenius law with zero activation energy and
temperature exponent and independent pressure law. Original mechanism admission
continues to own identities, elemental/charge and molecular-mass conservation,
and all species calorics. Ineligible mechanisms refuse at construction;
temperature-dependent or thermodynamic reverse rates are never frozen by this
selection. `ElementaryReactionV1::from_half_life` supplies eligible authored
half-lives without a second chemical rate convention.

Competing channels contribute their original hazards to each parent's total
removal rate. Original product stoichiometries determine yields without changing
multiplicity. The same positive owner prepares endpoint populations and
population-time integrals; channel events are `k_j * integral(N_parent dt)`.
This includes daughter ingrowth and preserves a small channel even when its
hazard no longer changes the rounded total hazard. An unrepresented positive
branch or event refuses. Counts are never inferred from a net inventory change.
An independent original-stoichiometry check uses the represented population
scale plus explicit counterproduct ULP allowances, rather than a relative
tolerance on large opposing event counts. This admits late daughter tails even
when both macroscopic counters round to the same value; the original caloric
and composition-reference-roundoff checks remain in force.
Separately, each interval charges one ULP of every channel counter, weighted
by the original source-combined reference reaction ΔU (or ΔH at fixed pressure)
and reference volume. The accumulated `first_order_counter_roundoff_j` reduces
the caloric root allowance, and that account plus the actual energy residual
must fit the original `energy_tolerance_j`. Large opposing cycle counts cannot
hide energetically significant quantization inside the reference-roundoff
diagnostic. This account is independent of arbitrary elemental formation
offsets; it is not a bound on all propagator arithmetic or source uncertainty.

The original reactor rebases directly onto the propagated populations and
independent channel extents, then performs its existing source-relative caloric
solve. Trace survivors are not reconstructed by subtracting nearly equal
amounts. Both adiabatic U,V and H,p entry points use this path, including explicit
energy initialization and restart. At fixed pressure, `V=nRT/p` and pressure
work remain the existing owner; the amount generator is independent of volume
because `V*k*(N/V)=kN`. Formation energy is included once. Final mass, element,
charge, energy and entropy gates are unchanged.

`maximum_step_s` bounds the interval between actual caloric and entropy checks;
`maximum_steps` bounds complete exponential intervals. Each interval has no
species time-discretization error for the declared constant generator. The
shared positive series, floating-point maps and source data still have their
own numerical/model limits; the method is not exact arithmetic or an accuracy
certificate. Initial/minimum steps and adaptive local-error controls belong to
the iterative methods. The balance reports zero Newton iterations and no local
time-error estimate. Thermal bath/finite-body exchange explicitly refuses this
method, and no continuous-interior source-support or entropy theorem is implied
by the checked states.

One immutable map is cached by exact interval duration and reused across
independent initial inventories. Shared node, edge, series, squaring and work
limits are explicit in `linear::FirstOrderLimits`; original cancellation and
no-partial-publication semantics remain. Checked map application allocates
bounded temporary arrays, unlike the unchanged allocation-free implicit
stages. The prepared `workspace_bytes()` report includes local binding vectors
but excludes the exponential owner's cached and temporary arrays.

Focused G1/G3/G4 regressions cover a half-life in one exponential interval,
equal-rate daughter ingrowth and calorimetry, rare competing channels, cyclic
events, trace survivors and sub-ULP parent losses, U,V restart and H,p work,
late daughter tails, excessive counter-energy quantization, ineligible laws,
numerical budgets, cancellation, zero time and deterministic
replay. They reuse the original source-backed NASA-9 reactor fixtures.

## Shared chemical sources for axial and planar conserved cells

`multicomponent::reactive::CellChemistry` binds the identical original
`MixtureGas` species bank and `IdealGasKineticsV1` mechanism, retaining the
reactor's selected time method and admitted source interval. Its axial
`[rho_i..., px, E_ref]` and planar `[rho_i..., px, py, E_ref]` updates share
the same conversion, original fixed-U,V reactor and endpoint reconstruction
as existing one-dimensional and network chemistry. Both planar kinetic terms
enter internal-energy admission; both momenta are preserved. Reaction changes
the stored reference gauge by the negative original species-reference change.
There is no additional heat-release source on physical chemical total energy.

The conservative endpoint is independently reconstructed at the original
reactor temperature. Caloric precision and the kinetic-energy consequence of
any molecular-mass residual must fit the declared energy tolerance. Each cell
and its accumulated report publish together only after the final cancellation
poll. Source, energy, work-counter, dimension and cancellation refusals leave
both unchanged. Numerical scratch is reusable after a refusal. Successful
advances allocate no heap storage.

The prepared reactor and cell source expose actual numerical heap capacities,
including original kinetics, constrained Newton, stage and amount buffers.
Immutable mechanisms, gas caloric preparation, fixed-size owners and bounded
stack scratch are explicitly outside that numerical budget. Higher-level
spatial owners can add these capacities to their own field workspaces before
reserving a second complete field.

Four focused regressions compare exothermic planar cells with the original
homogeneous reactor for backward Euler and SDIRK2, check the exact axial limit,
rotation and energy-datum invariance, and exercise late cancellation/replay,
source/path mismatch and complete cell/report rollback. This chemical-source
adapter alone supplies no spatial transport, molecular diffusion, ignition or
flame validation, and no global splitting-order or entropy certificate.

## Three-dimensional original-mixture Riemann and wall closures

`multicomponent::SpatialCell` and `SpatialPrimitive` extend the original
mixture owner to an orthonormal normal/two-tangent face frame. Canonical
species, all three momenta and full species-referenced total energy use the
same const-generic HLLC/HLL implementation as axial and planar consumers.
Star and HLL source admission includes both transverse kinetic contributions;
mixing their shear energy cannot escape the original caloric interval. The
original planar arithmetic and source-owned shock/rarefaction wall pressure
are retained. Reconstructed primitives retain source-computed temperature,
heat capacity and sound speed without another energy inversion.

Successful face calls use bounded stack scratch and caller output. Source and
precision refusals never clip temperature or renormalize species. Stationary
adiabatic slip walls have exactly zero species, tangential momentum and energy
flux and use the incident cell's original composition and calorics. Eight
focused tests cover bitwise planar reduction, material contacts with both
shears, tangential rotation/normal reversal/Galilean covariance, combined-shear
source refusal, variable-Cp and formation-gauge retention, original wall
pressure, precision/cancellation and the full 128-species bound.

This remains the Estimated conservative engineering lane. No moving
unequal-caloric interface pressure-equilibrium property, fully discrete
entropy/invariant-domain theorem, chemistry, spatial convergence or molecular
shock-thickness validation follows from a successful face evaluation.

### Conservative moving faces with evolving composition

`MixtureGas::moving_hllc_normal_flux`, `moving_hllc_spatial_flux` and their
reconstructed-primitive entries provide laboratory ALE transport through
translating internal faces. Every original species, all momentum components
and species-referenced total energy use the same face-relative HLLC/HLL
solution. The source reconstructs face-frame heat from each current
composition; both tangential kinetic terms participate in precision and
intermediate-state admission. No temporary frozen mixture or independent
passive-species flux replaces that owner.

The returned flux is `F - w_n U`. Its momentum and energy transformation uses
the actual ordered sum of species mass fluxes, preserving the original
per-species energy references and laboratory pressure work. A comoving
composition/shear contact transfers no species or tangential momentum, but
still transmits pressure traction and `p*w_n` work. Zero speed preserves the
stationary path; zero shear reduces to the existing axial moving-volume
kernel, which now shares this transformation. Successful calls allocate no
heap storage and retain the existing bounded species axis.

G0/G1/G3/G4 regressions exercise dimensional reductions, Reynolds transport,
comoving material contacts, vector Galilean covariance on HLLC/HLL, normal
reversal, source-reference invariance and caloric/cancellation refusals.
The spatial owner still supplies moving volumes, shared-face signs, discrete
geometric conservation and temporal integration. These fluxes do not add a
chemistry step or establish a pressure-equilibrium or entropy theorem.

## Complete methane chemistry with trace-preserving local extents

The existing prepared reactor now admits the complete 53-species, 325-channel
`kinetics::gri30::cantera_gri30` source, including its 16 originally irreversible
channels. `ContinuousAt300K` uses the source's explicitly declared upper-region
enthalpy/entropy continuity adjustment, retaining every original Cp law and
each species' own NASA seam. Its common admitted interval is 300--3000 K.
No reduced global methane reaction or added combustion heat source is used.

Reaction coordinates are rebased only after an accepted complete step.
Candidates reconstruct from that accepted composition plus local stoichiometric
extents; the original feed and fixed total internal energy remain untouched for
whole-call mass, elemental, charge, energy and entropy audits. This prevents a trace
intermediate from becoming the subtraction of large accumulated reaction-cycle
extents. It introduces no species clipping, positive trace floor, radical seed
or equilibrium projection, and rejected attempts retain their accepted origin.

The smaller Newton system remains algebraically equivalent to the full extent
equation. Exact nonzero stoichiometric rows and columns are prepared once,
without thresholds or rank truncation; compensated sums retain source order.
For GRI, the common factorization is 54 by 54 instead of 325 by 325. Its
stoichiometric matrix has 1227 nonzero entries out of 17225, reducing reduced
assembly product terms from 947700 to 83808, including the unchanged thermal
row. A focused release measurement of 40 identical source-derived equations
took 13.151 ms with the sparse traversal versus 90.728 ms with the dense
traversal; their 325 corrections were bit-identical. These are one-run
operation timings, not a whole-reactor or portable speed guarantee.

Every reduced correction must satisfy the original unscaled componentwise
linear backward-error gate. Up to four residual corrections are bounded.
If trace-channel cancellation still fails that gate, prepared storage permits
the existing scaled-partial-pivot full-extent solver, followed by the same gate
and at most two direct residual corrections. Physical and numerical acceptance
tolerances are unchanged. `IsochoricKineticsBalance::direct_fallbacks` reports
successful recovered Newton corrections, including those in later-rejected
time-step attempts. For GRI the fallback matrix/RHS/scales reserve 850200
bytes; the sparse column cache adds 22240 bytes beyond the sparse-row path.
All numerical buffers are prepared before advancing; no hot-path allocation
is introduced.

The independent development-only reproducer
`data/reference/thermochem/gri30_cantera32.py` verifies Cantera 3.2.0 and the
pinned source blob, then integrates the complete mechanism separately with
CVODES. Two native SDIRK2 regressions use CH4:O2:N2 = 1:4:15.04, one litre,
101325 Pa and no radical seeds. At 2 ms, the 1300 K induction case agrees with
the 1300.5483131534243 K reference within 0.00000181 K; the 1500 K ignition
case agrees with the 2600.310562321808 K reference within 0.0000561 K.
Both compare all 53 species within 0.5% plus 1e-12 mol, original physical
energy within 2e-8 J, pressure, independent atom/mass inventories and entropy.

The hot case uses 911 accepted steps, 154 rejected attempts, 11997 Newton
iterations and 1296 successful direct fallbacks at relative tolerance 2e-5.
Its observed mass residual is 1.11e-19 kg, maximum relative elemental residual
1.19e-15 and energy residual magnitude 2.85e-14 J. The induction case uses
110 accepted steps and no rejections. These are fixture-specific work and
accuracy observations, not certified error bounds. All 25 focused chemistry
tests pass, including original hydrogen, integration-order, conservation,
datum-invariance, cancellation and sparse/dense-equivalence regressions. The
complete focused native gas library also passes all 158 tests.

This extends homogeneous methane/natural-gas/nitrogen chemistry. It supplies
neither a fitted wood/PMMA pyrolysate nor soot, char oxidation, spray, real-gas
or experimentally validated flame/emission prediction. A finite-step entropy
check does not turn the source's irreversible channels into detailed-balanced
reactions, establish a fully discrete coupled entropy theorem, or give higher
order to surrounding pressure, transport, radiation or solid-source splitting.

## Shared original-mixture planar molecular faces

`multicomponent::dissipation::PlanarTransport` prepares the original
Maxwell--Stefan workspace and transport law on the exact evolving species
bank. Its shared private face kernel also serves the original axial and gas
network operators; there is no second species conversion, diffusion solver or
enthalpy convention. A planar face transports every diffusive species flux,
both momentum components and full referenced energy together. Its Newtonian
stress includes normal and transverse derivatives under the three-dimensional
trace convention with zero out-of-plane velocity. Energy includes partial
referenced species enthalpy, Fourier heat and both components of stress work.
No separate chemical or viscous heat source duplicates those transfers.

Conservative inputs use the original caloric inversion. Already reconstructed
primitive inputs use the original source admission without another inversion.
Stationary adiabatic slip walls retain normal viscous traction and exactly zero
species, shear and energy flux. Output is scratch until success; callers own
spatial gradients, conservative divergence, requested-duration transaction and
the complete boundary impulse. The reported diffusion scale is a conservative
engineering step-size estimate, not a spectral or invariant-domain certificate.

Prepared numerical capacity includes the binary matrix and all original
Maxwell--Stefan solve buffers. Successful kernels allocate no new scratch.
Immutable gas/species storage, bounded stack scratch and the caller's transport
law storage/evaluation allocations are outside that budget. Material-card
providers retain their original source-query allocation and validity contract.

Five focused regressions cover independent binary diffusion with unequal Cp,
Fourier/full-tensor/enthalpy flux, primitive/conservative agreement, tangential
frame and formation-gauge covariance, exact axial reduction, slip-wall loads,
exact-zero species, source/axis refusal and final cancellation/replay. Original
axial/network regressions remain authoritative for their geometry and balances.
Soret/Dufour coupling uses the explicit source-law extension above. These face
fixtures establish no nonideal transport, spatial accuracy or combustion
validation. Explicit held
reservoir faces use the source-owned boundary energy convention below.

### Shared original-mixture molecular faces in three dimensions

`multicomponent::dissipation::SpatialTransport` extends the same prepared
source and const-generic molecular kernel to all three momentum components.
The original axial and planar owners retain their public interfaces and use
that same kernel with zero or one transverse component; no diffusion solver,
species conversion, caloric reference or constitutive law is duplicated.
Three-dimensional conservative inputs remove all kinetic components through
the original source inversion. Primitive inputs use the original complete
three-momentum admission without repeating the temperature solve.
Every face receipt retains `dynamic_viscosity_pa_s` from that original law
evaluation. A composing SGS owner can check its eddy/molecular ratio against
the actual face coefficient without querying a stateful source a second time.

Each face supplies its positive normal separation and both tangential rows of
the velocity gradient in the original orthonormal face frame. The Newtonian
traction includes both transverse dilatations in the normal component and
both normal/tangential derivatives in each shear component. All three stress
work terms accompany original referenced species enthalpy diffusion and
Fourier heat in one total-energy flux. Stationary impermeable adiabatic slip
walls retain both transverse dilatations in the normal viscous load while
returning exactly zero species, shear and energy transfer. A field owner must
still supply physical gradients, divergence, time integration and complete
boundary exchanges.

The prepared binary matrix and Maxwell--Stefan numerical heap capacities are
unchanged by the extra momentum direction. Successful kernels allocate no
additional scratch; caller-supplied law allocations keep their own contract.
The reported engineering diffusivity estimate includes the complete
three-dimensional Newtonian row couplings and must be combined with all
spatial inverse lengths by the field owner. It is not a spectral,
invariant-domain or fully discrete entropy certificate.

Six focused regressions independently reconstruct binary species/Fourier/full
stress-work fluxes, compare conservative and primitive inputs, rotate both
tangent axes and the complete derivative tensor, check Galilean and formation
reference covariance, retain the planar reduction bit for bit, and verify
slip-wall loads, a resolved `1e-180` trace species, an exact absent species,
source/axis/gradient/caloric refusal and final cancellation/replay. They do
not establish spatial convergence, physically resolved shock/flame thickness,
Soret/Dufour effects, no-slip walls or experimental combustion validity.

### Original-source energy effort at physical reservoir faces

Both planar and spatial face owners expose `reservoir_face_from_primitives`
for an explicitly selected exterior trace at a physical boundary. The original
source admits both traces. Maxwell--Stefan, traction and Fourier transport
retain the shared two-trace midpoint approximation, while species energy uses
the bath's referenced enthalpies and stress work uses the bath's velocity.
A stationary bath therefore receives zero mechanical work even with nonzero
viscous traction. The supplied distance is cell-center to physical boundary,
not the spacing between two cells. The existing three-dimensional gradient
layout and all interior-face source conventions remain unchanged.

The kernel does not evolve the exterior reservoir. The field owner must retain
that state, its geometry and signed boundary transfers consistently with its
advective boundary condition. A focused regression checks both exterior-side
orientations, original unequal-Cp bath enthalpy with large formation datums,
all three work components, unchanged species/traction fluxes, and zero work at
a stationary bath in both dimensions. Spatial accuracy and finite-reservoir
storage require their respective field and coupling owners.

### Three-momentum chemistry on the same prepared cell source

`CellChemistry::advance_spatial` extends the same prepared axial/planar owner
to canonical species, px/py/pz and referenced total energy. All three momentum
components remain fixed, and all three kinetic terms enter the original-source
internal-energy target and independent endpoint reconstruction. The original
formation-reference change is applied once, with no extra heat-release source
on chemical total energy. Dimension-specific public entry points share one
const-generic conversion and acceptance kernel.

The existing atomic cell/report publication, exact source identity, caloric
interval and bounded-work requirements remain unchanged. Actual prepared
reactor capacities are available to the enclosing spatial memory budget.
Focused dimension-reduction/rotation and late-cancellation/replay checks cover
complete chemical energy, all momenta and reusable prepared storage. This is
the local fixed-volume chemical source; it does not itself establish a spatial
entropy theorem, resolved flame structure or experimental combustion validity.

### Compatible composition through density and compression waves

`multicomponent::reconstruction::mass_fraction_slopes` supplies the same
bounded, allocation-free composition limiter to axial and Cartesian mixture
consumers. Traces reconstruct density and every canonical mass fraction
separately, then recover partial densities as their product. Conservative
cell averages and their original species equations are unchanged. A single
MC factor scales the entire centered fraction-slope vector; constant inert
fractions and other linear composition constraints survive density variation
to floating-point roundoff. No dependent remainder species, normalization,
trace seed or independent species correction is introduced.

The relative-only 32-epsilon roundoff guard prevents division noise in an
otherwise constant fraction from flattening all composition slopes, while
retaining resolved gradients even for a 1e-200 trace constituent and preserving
exact absent species. The axial owner uses actual nonuniform cell-centre
distances. Its source-domain contraction scales density, every fraction,
velocity and pressure together before the original thermodynamic admission.
The additional density reconstruction scalar is prepared once per cell and
does not add a second conserved density equation or per-step allocation.

Focused compression-wave tests on uniform and nonuniform 80-cell meshes
retain a uniform inert fraction to 2.23e-16 and an independent weighted
composition constraint to 4.45e-16 while closing every species and energy
balance. These verify compatible reconstruction; they do not establish a
fully discrete composition-maximum principle or entropy theorem.


## Quasistatic moving compact gas volumes

\`multicomponent::network::moving\` adds the first explicit moving-boundary
thermodynamic owner for species-bearing compact gas storage. A
\`Plenum::quasistatic_adiabatic_to\` candidate preserves every canonical species
mass while changing physical volume. It solves the original fixed-composition
NASA-9 ideal-gas isentrope, \`Delta s(T,p)=0\` with \`p=rho R_mix T\`, through
the same borrowed caloric source used by the wall/characteristic owners.
Safeguarded Newton/bisection is bounded by the existing caloric-domain work and
temperature limits; unsupported compression/expansion refuses rather than
clipping temperature or pressure.

The receipt reports source temperature/pressure endpoints, mechanical boundary
work on the gas, physical internal-energy change, extensive entropy defect,
per-species mass residual and the referenced-versus-absolute energy-gauge
residual. Because the closed species inventories are unchanged, formation
references cancel from work exactly in the continuum. Boundary work is not
reclassified as heat and no \`p*Delta V\` endpoint approximation is used.

`MixtureNetwork1d::quasistatic_adiabatic_plenum_to` accepts that same endpoint
directly into a running connected chamber between transport steps. Every branch
field, other plenum, connection and clock stays unchanged. The next transport
step reads the accepted chamber volume and source state for both port fluxes
and its acoustic CFL restriction. A single branch-field pass admits old/new
combined inventories before publication, including refusal if total network
storage would overflow. The returned `NetworkVolumeReport` separates mechanical
work from transport/heat and reports the combined species and energy balance.
An actual accepted volume change invalidates the obsolete flow/port reports;
an exact no-op, refusal or cancellation retains them. No network reconstruction
or source-gauge conversion is required.

These operations are deliberately clock-free and quasistatic. Network use
holds ports closed during the volume change and leaves the schedule and
accuracy of composing that change with transport to its caller. They do not
implement a finite-speed moving-wall Riemann solver. Piston inertia, prescribed
motion, shock generation during boundary motion, heat leakage and moving
geometry/GCL remain explicit separate owners. Three focused endpoint tests cover an independent
constant-Cp isentrope and reversible round trip, variable calorics with +/-1e18
formation-energy gauges, and source-domain/cancellation atomicity. Three
additional network regressions cover compression-driven branch flow with an
independently integrated physical work balance, reversible volume cycles, and
state/clock/report preservation under early or late cancellation, source
refusal and exact no-op.

## Finite-time compression coupled to network transport

`MixtureNetwork1d::step_with_volume_rates` evolves prescribed per-plenum dV/dt
simultaneously with branch and characteristic-port flow through the original
SSP RK2 stages. Species inventory and referenced internal energy advance
extensively, with source-owned pressure work `-p dV/dt` accounted separately
from transported energy. Accepted endpoint volumes and pressure work are
available per plenum, and `NetworkStepReport::boundary_work_on_gas_j` closes
the combined network energy balance with that same pressure quadrature.

Positive representable volumes and acoustic/volume-rate CFL bounds gate atomic
publication. All-zero rates preserve the rigid-step arithmetic; failure retains
fields, volumes, clocks and accepted receipts. Work/report buffers are allocated
at construction and reused. G1 tests compare 20/40/80-step compression and
expansion ladders with independent constant/linear-Cp isentropes. Other focused
tests check simultaneous port enthalpy/species/work closure, exact rigid
identity, geometry/source-domain refusal and early/intermediate/final cancellation.

The compact perfectly mixed chamber approximation remains. Rates are constant
over each accepted step. Spatial piston shocks, branch-wall motion and piston
inertia require other physical owners; this is not a spatial ALE solver.

### Closed-network entropy admission with moving chambers

`MixtureNetwork1d::step_with_entropy` and
`step_with_volume_rates_and_entropy` apply one original-source entropy account
to every branch cell and finite plenum inside the closed network. The source's
temperature-dependent `Cv/T` integral and ideal-mixture composition entropy
use actual local states and physical volumes. Internal characteristic-port
exchanges are already inside this boundary and add no exterior entropy supply.
Prescribed plenum pressure work retains the native energy quadrature and
carries no entropy. Unconnected branch ends must be stationary reflecting
walls; isolated periodic branches are also admitted.

The optional gate uses the existing bounded SSP RK2 retry/publication loop.
An entropy-decreasing candidate retries all owners from the accepted state;
source-domain failures, exhaustion and cancellation publish no fields, volumes,
clock or accepted work/exchange reports. The returned duration is the actual
accepted duration. Every species inventory defect is independently bounded,
with only its numerical source-reference entropy contribution removed from
the gate and reported separately. No energy, species or temperature correction
is applied. Original and final caloric states share one source/reference basis;
moving-volume increments retain extensive mass differences and stable density
remainders rather than subtracting large absolute entropies.

Two bounded temperature arrays are prepared once per entropy-enabled call;
retries reuse them and native numerical storage. Ordinary step arithmetic and
all-zero-volume-rate identity are retained. G0/G1/G3/G4 regressions compare
polynomial-Cp network entropy with independent extensive inventories, moving
mixture chambers with reversible isentropes and pressure work, large independent
source entropy/formation datums, whole-network retry, exhaustion and final
cancellation. This is measured finite-step admission of the native engineering
flux. It does not establish an entropy-stable semidiscrete port theorem, and a
negative semidiscrete rate can exhaust retries. The tolerance is per accepted
call; accumulated multi-step allowances remain caller-owned. Open reservoirs,
chemistry, thermal sources and moving branch walls need their own entropy inputs.

### Prescribed reservoir entropy through complete moving networks

`MixtureNetwork1d::step_with_reservoirs_and_entropy` and
`step_with_reservoirs_and_volume_rates_and_entropy` extend that same native
transaction to prescribed stationary pressure reservoirs. One global source
basis covers every stored branch cell and finite plenum. The existing
`ExteriorEntropyEffort` contracts the actual signed SSP species, momentum and
referenced-energy transfers with the original reservoir entropy variables and
advective entropy potential. For stationary reservoirs this reduces to their
heat/chemical effort on transported total enthalpy. Pressure-flow work is
already included in enthalpy; moving-plenum `-p*dV` work remains a separate
mechanical energy exchange with zero entropy supply.

Topology selects external physical ends. Internal finite-plenum ports and
periodic faces contribute no boundary entropy or external species throughput.
Unconnected ends without reservoirs must be reflecting walls and carry zero
species/energy at both stages. Nonzero transfer of a species with zero reservoir
activity refuses at either stage; opposite transfers cannot hide it. Cached
efforts are prepared once per call in the complete network's fixed reference
basis. Retries reuse them and the native numerical workspace.

`ExteriorMixtureEntropyReport::boundary_entropy_input_j_k` sums all left and all
right physical ends, respectively, in retained branch order. One global gate
checks state entropy minus this supply; it does not sum branch-local allowances.
Each species continuity check subtracts the actual boundary transfer and scales
by original inventory plus absolute throughput from both stages and every
external end. Only the numerical continuity reference contribution is removed;
physical entropy datums carried by entering/leaving species remain in state and
boundary reports. Receiving reservoirs retain their own thermodynamics, so the
reported generation includes irreversibility of the imposed reservoir coupling.

Focused regressions independently invert polynomial-Cp states, reconstruct
reservoir chemical/enthalpy entropy from accepted native transfers, close moving
work, admit entropy-decreasing outflow with positive combined generation, retain
large incoming-species entropy datums, and check whole-network retry, rollback,
late cancellation, topology, zero activity and closed/zero-motion identities.
The accepted duration may be shortened. No field, geometry, clock or accepted
receipt publishes on failure. This remains measured finite-step admission of
the existing characteristic engineering flux, without a fully discrete entropy
or invariant-domain theorem or finite-reservoir depletion model.

### Frozen-composition chambers with wall and gravity coupling

`dynamics::network::GasNetwork1d::step_with_volume_rates` supplies the same
prescribed compression/expansion to the original frozen-composition owner.
Its SSP RK2 stages evolve extensive mass and source-referenced internal energy
on the changing physical volume, with `-p dV/dt` supplied by the original
variable-Cp caloric source. `accepted_plenum_boundary_work` uses the same
`PlenumBoundaryWork` receipt as the evolving-mixture network. Construction
allocates the stage/report buffers once; all-zero rates preserve rigid-step
arithmetic. Positive represented volumes and stage acoustic/volume-change CFL
limits precede atomic state, volume, clock and report publication.

`step_with_walls_and_volume_rates` composes the existing implicit finite-wall,
shared-wall and branch-loss sources at the actual accepted endpoint volumes.
Wall geometry bindings update only on complete acceptance; retries rebuild
thermal candidates from the accepted wall inventory and shortened geometry.
Heat/work stay separate in gas, finite-wall and gravity balances. Compact
gravitational storage uses the candidate volume and unchanged single-potential
location, so compression cannot create gravitational mass. Authored wall
conductance, area and thermal inventory are not inferred from a volume rate.
Transport/thermal splitting remains first order despite second-order transport.

Focused G0/G1/G3/G4 tests cover constant and linear-Cp compression/expansion
against independent entropy integrals, connected pressure-wave discharge and
work conservation, analytic finite-wall/shared-wall endpoints, gravity-datum
invariance, formation-gauge and rigid-path identity, volume CFL, and early/late
cancellation and source refusal without partial publication. These numerical
conformance tests do not validate real engine geometry or a fully discrete
entropy theorem.

### Moving chambers with reacting source operators

`ReactingNetwork1d` composes prescribed plenum volume rates with finite-rate
chemistry, molecular species/heat transport, viscosity and finite-wall storage.
The existing split is W/2 -> V/2 -> D/2 -> C/2 -> moving transport -> C/2
-> D/2 -> V/2 -> W/2, omitting absent operators. Source halves before transport use the
original volume; halves afterward use the actual accepted endpoint volume.
Boundary work has its own report field and is subtracted exactly once from
referenced and physical-energy residuals.

Prepared wall, molecular and valve operators keep accepted geometry separate
from unpublished interval volumes. Complete-step acceptance publishes the new
bindings, while shortened-duration retries restore original fields AND volumes.
Final cancellation preserves accepted geometry, reports, valve exchanges and
wall storage. All-zero rates retain the original rigid path, including valves.

Valve halves use the actual original/final chamber volumes in the unchanged
source-owned nozzle and finite donor-blowdown solve. Each half holds its stage
volumes fixed and transfers identical species and finite enthalpy to the
receiver; only moving transport supplies mechanical pressure work. Source
receipts sum the two valve intervals' storage defects, excluding intervening
compression, chemistry and heat. Public standalone valve calls still require
the accepted geometry, and failed trials cannot rebind it.

Authored wall conductances, contact areas and compact mixing lengths remain
fixed inputs; a volume rate does not determine their changing geometry.
The compact mixing and existing reaction/source-splitting approximations remain.
Focused regressions exercise live reaction/transport/compression/wall energy,
duration restart, final cancellation and reuse, accepted geometry rebinding and
exact zero-motion identity. Moving-valve cases additionally reconstruct each
tank's species and energy change from finite enthalpy exchanges and pressure
work, close physical gas-plus-wall storage with all sources active, and check
complete-duration retries and cancellation/restart. These new valve-composition
regressions are added but have not been executed in this resource-constrained
session. No new global temporal-order or entropy-stability claim is made.

## Original material calorics in finite gas walls

`MaterialWallCaloric` binds the existing `PreparedHeatStorage` to finite
branch-cell and plenum walls. Its explicit caloric constraint declares
fixed-volume internal-energy capacity or the stationary-solid
negligible-expansion approximation; Cp is never silently converted to Cv.
Original source lines, fixed constituent amounts and temperature support remain
attached. Source integrals determine actual storage; a temperature contact's
capacity is only its refreshed tangent cache.

The same binding's `from_thermal_law` retains an original
`fs-thermal::ThermalLaw`, its physical reference cell, fixed mass and explicit
supported energy interval. Both owned phase charts and adapters borrowing
prepared material sources can be bound without copying a fitted heat capacity.
`NetworkWallOperator::with_material_heat_storage` supplies either source to the existing
ordinary, molecular, reacting and valve-composed network transactions.
Developed pipe-wall transport continues to query current gas composition,
viscosity and conductivity; it replaces only drag/conductance, preserving the
original wall source and evolving stored energy.

The source step solves one bracketed backward-Euler equation in transferred heat:
`Q = dt G [T_wall(e0-Q/m)-T_gas(E0+Q/V)]`. The original gas calorics recover gas
temperature and the original wall law supplies T(h/u), including exact latent
plateaus. Monotonicity belongs to the fixed-conductance case; a source-dependent
endpoint conductance makes no uniqueness claim. The same safeguarded scalar owner uses exact source tangents when
available and bisection otherwise. No apparent capacity or separate latent-heat
source is added. Actual
represented gas energy, wall storage and heat-flux residuals must all meet the
existing joule allowance. Source/geometry/range/work/cancellation failure returns
no candidate; the original complete-network transaction publishes every gas and
wall state together. Source storage is immutable and numerical work is reused.

The law's storage basis remains explicit: internal energy describes an admitted
fixed-volume source; enthalpy retains the material's own pressure-work convention.
The combined report sums the selected thermal coordinates and does not relabel
wall enthalpy as rigid-vessel internal energy. Reference heat-transfer geometry,
uniform-wall applicability, imposed conductance and any external pressure work
remain caller/source responsibilities. Melting-wall motion, loss of mechanical
support, wall sharing between cells, nucleation, and spatial wall conduction are
not inferred by this finite-store coupling.

Focused tests check independent nonlinear calorimetric references, original
source-only substitutions, reaction/domain/cancellation rollback, latent storage
and full plateau traversal in both heat directions, accepted-history rollback
and replay, continuing molecular/wave energy conservation, and the original
borrowed solid material chart in a reacting network with current gas-card losses.
Native execution remains unverified in the authoring worker: the unchanged
asupersync dependency was SIGKILLed during code generation under its memory cap.

Four further G1/G3/G4 cases exercise implicit material conductivity against
independent latent-plateau quadratic roots in both heat directions, a valid root
beyond the frozen-flux bound, variable-capacity endpoint energy, source-only
changes, constant-conductance and disabled-film limits, continuation and
source/work/cancellation refusal. The reacting-network regression runs both time
evaluations. These new Rust cases remain unexecuted; six changed Rust files
passed pinned syntax parsing and scoped whitespace checks.


## Pressure-actuated network valves

`Orifice::with_pressure_opening(PressureOpening)` adds a reversible quasi-static
differential-pressure aperture to the existing finite-plenum valve owner.
Finite typed thresholds satisfy `0 <= cracking < full_open`: aperture is zero
through cracking, linear between thresholds, and capped at one. It multiplies
the declared maximum `Cd * opening * area` exactly once. Existing
`ValveDirection` enforces exact reverse closure.

Initial and implicit endpoint evaluations use current source-decoded pressures.
The actual nozzle pressure ratio, source-owned variable calorics and choking,
adiabatic donor blowdown, and equal-and-opposite species and energy transfers
remain authoritative. Cracking pressure is never subtracted from nozzle
backpressure. Closed endpoint trials bracket the transfer instead of accepting
threshold overshoot under a loose absolute mass tolerance. Invalid thresholds,
nonfinite differences, and unresolved nonzero aperture or area refuse.

This is an Estimated algebraic aperture law without hysteresis, poppet inertia,
spring-energy storage, leakage, calibration, or a certified safety-relief claim.
Fixed-opening construction retains its previous behavior. Focused
`pressure_actuated` regressions exercise analytical choking/un-choking,
thresholds and directions, finite-tank blowdown and energy/species balance,
source rollback, and pressure-wave propagation in the original gas network.
Their presence is not an executed validation claim.


## Prescribed species-bearing pressure reservoirs

`MixtureEuler1d::step_with_reservoirs` accepts stationary species-bearing end
cells in the same source species order and per-species energy-reference gauges.
It recomputes the existing characteristic fluxes from the actual interior cell
at both SSP RK stages and every halved retry. Reservoir inflow carries the
source composition and stagnation enthalpy; outflow carries the evolving
interior composition and total enthalpy. Reservoir pressure supplies inflow
stagnation pressure or outflow backpressure through the original variable-
caloric shock/fan and choking owner.

An absent endpoint preserves its declared boundary; two absent reservoirs
execute the original step exactly. Periodic intervals with reservoirs, moving
reservoirs, invalid dimensions and source-domain violations refuse before cell
or clock publication. Borrowed `Cell` values carry no owner identity token:
callers retain responsibility for matching canonical species order and energy
gauges, as with `MixtureGas::reservoir_flux`. Stack endpoint buffers add no
per-step heap workspace. Prescribed states are fixed during one step and are
not depleted; finite-reservoir feedback remains with the existing network owner.

Six focused `reservoir_boundary_tests` regressions exercise independently
reconstructed predictor/corrector fluxes, variable-Cp donor enthalpy, pressure/
temperature/composition drives, end reversal, formation-reference invariance,
CFL retry replay, exact absent-reservoir parity and cancellation/refusal rollback.
Their presence does not establish executed native validation.


## Prescribed mixture pressure reservoirs with finite-volume motion

`MixtureNetwork1d::step_with_reservoirs` accepts one optional stationary left/right source pair per branch, in canonical species order and energy reference gauges. Raw `Cell` source compatibility remains caller-owned; species dimensions, stationarity and source-state validity are admitted. Only unconnected physical ends may have prescribed reservoirs. Characteristic fluxes refresh at both SSP RK2 stages and every retry, preserving reservoir inflow species and stagnation enthalpy, and interior outflow species and total enthalpy. External fluxes enter the global boundary balance; finite-plenum exchanges remain exact opposites. An all-None source slice keeps the original step. Fields, clocks and accepted reports publish atomically. Construction-time primitive caches and existing endpoint workspaces avoid per-step numerical heap allocation.

`step_with_reservoirs_and_volume_rates` applies prescribed exterior sources and finite-plenum motion in the same native SSP RK2 transaction. Each stage uses the evolving chamber volume, composition and energy alongside fresh branch-end characteristic states; mechanical `-p*dV/dt` updates extensive energy beside conservative port transport. Reports separate exterior enthalpy transfer from accepted plenum boundary work. All-None sources preserve the volume-rate lane; all-zero rates preserve rigid-reservoir arithmetic, and every retry replays both drives before one shared publication point. Reservoir controls are not yet propagated through reactive or wall-source composed network wrappers.


## Pressure reservoirs through coupled reacting and molecular network owners

The prescribed species-bearing pressure reservoirs of `MixtureNetwork1d` now
propagate through the complete reacting and molecular network transactions
instead of requiring a caller-side split. `ReactingNetwork1d::step_with_reservoirs`
and `step_with_reservoirs_and_sources` compose the original C/A/C chemistry
with reservoir-driven characteristic transport and optional wall/valve halves.
`step_with_reservoirs_and_volume_rates` evaluates the same exterior pressure
drive simultaneously with finite-plenum pressure work. The molecular owner
likewise exposes `DissipatingNetwork1d::step_with_reservoirs` and
`step_with_reservoirs_and_sources`, retaining D/2-A-D/2 and optional source
halves in one complete-duration retry.

Reservoir states remain prescribed and undepleted during one call. Every
transport retry recomputes their characteristic flux from the actual branch
stage state. Advected reservoir mass/enthalpy appears once in the transport
boundary report. When molecular evolution is selected, its two halves also
use those baths and retain their distinct diffusive species/heat transfers.
Chemistry reference changes, wall heat, valve-internal exchange and
moving-volume work keep their separate accounts. A shortened acoustic step
restarts all chemistry, molecular,
wall and valve owners from accepted storage before replaying the pressure drive.
Cancellation or source refusal publishes none of those states or receipts.

Focused regressions drive a new species from a high-pressure reservoir through
the reacting network, including simultaneous plenum compression, and through
the molecular D/A/D owner. They reconstruct species and physical-energy changes
from the complete reports. These tests are authored but native execution remains
subject to the repository's existing build-resource limits.

## Gravity-aware finite-plenum networks

`GasNetwork1d::new_with_gravity` couples the original isothermal and isentropic
branches through one explicit `NetworkGravity` convention. Each branch retains
its signed local acceleration and reference position, with
`phi(x) = branch_offset - g*(x-x_ref)`. A branch without gravity has constant
potential. Every rigid compact plenum has one declared specific potential;
all incident physical faces must meet it. Admission subtracts branch/plenum
datums before adding the local head, with tolerance based on that physical
difference rather than a large arbitrary absolute zero. Missing declarations,
mismatched heights, non-finite or unrepresentable storage refuse.

Connected characteristic ports reconstruct the physical hydrostatic face from
the original retained reference and current cell perturbation. Volume-averaged
end-cell pressure is not substituted. Identical admitted resting traces reuse
the original source evaluation. `Euler1d::hydrostatic_reference_face` exposes
the retained face state for consistent compact-volume initialization. Each
port still evaluates one source-owned characteristic flux transferred with
opposite signs. Gravity work uses the original branch face-mass-flux identity;
no extra junction work term is introduced.

Gas reports retain branch kinetic/internal energy and compact internal energy
and now explicitly account for branch gravity work. `NetworkGravityReport`
separately balances gas plus midpoint branch/compact potential storage against
external gas energy and `phi*mass` transport. All potential terms are accumulated
after subtracting the first branch offset. The exposed reference and centered
storage/transport preserve resolved work under a large common datum; the
absolute-potential residual is `R_centered + phi_ref*R_mass`. Absolute source
formation energy additionally uses the existing `e_ref*mass` convention. Direct
subtraction of huge absolute storage levels is not used for the conditioned
residual. No per-cell arrays are allocated during stepping.

Six native regressions cover exact joined isothermal/isentropic rest with both
reconstruction orders, transient exchange among three elevations with independent
storage reconstruction, external potential-energy transport, bitwise flow and
conditioned-balance invariance under a representable `2^40 J/kg` potential shift
and `1e18` formation datum, physical height/refusal checks, and rollback at the
final publication poll. All six new tests and all 80 native dynamics tests pass.

A compact node cannot represent a hydrostatic column spanning different
elevations. Original constant-area/constant-acceleration frozen-ideal-mixture
assumptions remain. This adds neither arbitrary entropy stratification, gravity
for reacting or diffusive networks, moving nodes, self-gravity, nor an
entropy/positivity certificate.

## Coupled atomic ionization and collisional thermal evolution

`chemistry::ionization::advance_atomic_ionization_with_collisions` evolves the
original ground-state `e + A <=> A+ + 2e` channel and electron/heavy thermal
exchange simultaneously in the existing rigid-parcel owner. It consumes
`energy::ElectronHeavyCollisionSources`: original neutral momentum cross
sections, optional Coulomb model, and explicit charged-fraction, massive-target
and Maxwellian quadrature controls. The original Coulomb receipt must match
the initial population. Later evaluations explicitly refresh the same source
law at the current amounts, density, and separate temperatures. Recombination
requires the neutral source even when that neutral is initially absent.

For reaction extent xi and cumulative electron-to-heavy transfer Q, the
electron budget is `Cv_e*(ne0+xi)*Te=Cv_e*ne0*Te0+Pe*t-chi*xi-Q`, with
`Cv_e=3R/2`. The heavy source-difference caloric budget is `H(Th)=Ph*t+Q`.
Thus chemical storage, changing electron heat capacity, external work and
internal heat transfer enter once. Inert species retain their original heavy
calorics. Exactly absent electrons remain absent and cannot receive power.

Each backward-Euler substep refreshes actual neutral and Coulomb conductances
at its start, then solves a scalar reaction residual with an inner monotone
thermal solve. The analytic constrained derivative includes both changing
electron population and heavy-pool thermal response. Safeguarded brackets
admit supported interior states without clipping temperatures or abundances.
Full/two-half-step estimates control this first-order coefficient lag and the
combined reaction/thermal evolution. Both temperatures, all reaction amounts,
and separate neutral/Coulomb transfers enter the error norm. Source admission,
physical mixture entropy and independent pool-energy balances gate publication.
Nonlinear trials allocate no state or species workspace; source snapshots are
refreshed between substeps. Cancellation leaves caller-owned input unchanged.

G1/G3/G4 regressions compare three refinement levels with an independent
NIST-BEB/Landau/analytic-neutral ODE calculation, exercise signed exchange and
held powers under formation-energy gauge changes, create an initially absent
neutral through recombination, and refuse stale source bindings or cancellation.
The native EM integration carries cell-volume heating through this same joint
update and refreshes transport from actual changed populations and temperatures.
These are numerical checks of the declared classical Maxwellian ground-state
model. They do not establish experimental discharge validation, radiation or
photoionization, excited-state kinetics, a continuously fixed-voltage transient,
spatial charged transport, or a general multireaction plasma network.

## Finite-plenum thermal equilibrium in gravity networks

`GasNetwork1d::step_with_thermal` couples the original frozen-gas transport to
a prepared `PlenumThermalOperator`. Each compact volume has one explicit
`PlenumThermalContact`: adiabatic, a prescribed-temperature reservoir,
authored finite capacity, or retained original `MaterialWallCaloric` storage.
The latter keeps its nonlinear h/u integral, original identity/geometry and
latent inventory; a plateau is not replaced by apparent heat capacity.
Distinct contacts represent distinct physical wall masses.

The existing implicit gas/wall source now borrows either the original frozen
or species-carrying gas calorics. Both adapters retain their actual source,
energy gauge and inversion; no equivalent species bank or constant-gamma
substitute is constructed. Existing mixture entry points and source-dependent
wall-conductance behavior remain on the same equation owner. Exact
source-representable thermal fixed points avoid introducing inversion
tolerance into zero heat flux. Passive source brackets, original source
domains, actual represented heat, energy residuals and bounded work remain
authoritative.

The wave owner first prepares its ordinary SSPRK2 candidate. Backward-Euler
plenum heat then uses that trial's actual accepted duration, including CFL and
source retries. This is first-order thermal/transport splitting. All gas
states, finite-wall states, clock, port exchanges and reports publish together
only after the final poll. Failed later contacts and final cancellation retain
both accepted owners; prepared source buffers are reused. Amount and volume
are fixed during heat exchange, so it adds neither momentum nor gravitational
work.

Gas balances subtract heat exactly once. `NetworkThermalReport` separately
records the original wall source and closes gas plus finite-wall storage
against external gas/reservoir transfer; with gravity it also includes
centered potential storage and external potential transport. The complete
absolute-potential residual is reconstructed from the complete centered
residual plus `phi_ref*mass_residual`. Formation energy retains the original
separate convention.

Seven focused native regressions cover exact stiff hydrostatic/thermal
equilibrium, independent variable-Cv finite-wall roots and long-time
equilibration, heating/cooling without temperature crossing, actual accepted
duration, nonlinear original wall storage, latent inventory, pressure-driven
gravity flow with independently reconstructed energy, large potential/formation
datum invariance, final cancellation/reuse and later-source rollback.
The seven new cases, all 87 dynamics tests and 14 existing shared-wall
regressions passed against the original production crate. The material-enabled
production path compiled. Its supplemental gate initially rejected three shared
fixtures with undeclared T/phi quantity types; their validity declarations now
match the original source support and strict query types. The corrected
material-only cases remain unverified: the rerun was killed by the shared-memory
limit during test compilation.

This remains a compact-volume, frozen ideal-mixture engineering model with
explicitly authored conductance. It supplies neither a validated heat-transfer
correlation, resolved wall conduction or shared wall mass, moving-volume work,
gas phase change, nor a second-order or entropy-stability certificate.

## Distributed wall equilibrium in gravity-aware gas networks

`NetworkWallOperator` extends the same prepared wall owner to explicitly
selected finite-volume branch cells and compact plenums. The original
`PlenumThermalOperator` name is an alias of this owner; existing plenum calls
retain their behavior. `BranchWallContact` binds original linear/Darcy drag
and reservoir, finite sensible, nonlinear or latent thermal storage. Original
source, complete mesh geometry, cell volume and distinct wall amount remain
explicit. A shared physical wall requires its own thermal owner rather than
duplicated stored energy.

`GasNetwork1d::step_with_walls` uses the existing transactional transport
step. The shared caloric wall kernel acts on unpublished branch and plenum
candidates over the actual accepted retry duration. Transport/thermal coupling
is first-order splitting. Drag converts kinetic to internal energy while
retaining total gas energy; only accepted wall heat changes that total.
Fixed gas amount and cell volume mean no additional gravitational work during
the wall source. Original variable heat capacity and latent h/u inventory are
retained without constructing an equivalent gas or apparent wall heat capacity.

Each accepted branch `StepReport` exposes signed axial wall impulse and
wall heat and includes both in its residual. Impulses along differently
oriented branches are not summed into a fictitious scalar network momentum.
The shared source ledger closes gas plus finite-wall storage and, where
declared, centered gravitational potential. Absolute-potential residuals
retain the existing reference-times-mass-residual conversion. Numerical
candidate arrays are allocated at preparation and reused. Gas, walls, clocks,
port exchanges and reports publish together after the final cancellation
poll; a later source refusal retains all accepted state.

Seven native regressions cover exact hydrostatic/thermal rest, analytic linear
and Darcy drag, independent variable-Cv finite-wall exchange, heat-driven
gravity flow with potential/formation-reference invariance, retained latent
inventory, final cancellation and later-source failure. The complete dynamics
suite passed 94 tests, including those seven; all 14 original shared-wall
regressions also passed (108 distinct tests). These used the original
production crate with `experimental-euler`, without the optional
`material` feature.

This remains a frozen ideal-mixture passage/plenum model with explicit wall
laws. It adds no gas phase change, resolved wall conduction, shared thermal
wall network, moving-volume work, calibrated heat-transfer correlation,
second-order splitting or entropy-stability certificate.

## Shared finite walls across frozen gas networks

With `experimental-euler`, `dynamics::network::NetworkWallOperator::with_shared_walls`
binds one finite original solid inventory to several branch or plenum gas
control volumes. A `SharedWall` owns that storage once; its `SharedWallPort`
values identify the original gas cells and their complete nonnegative
conductances in W/K. Sharing heat never mixes their mass or species. Local
contacts must be adiabatic without competing finite storage, and a gas target
may belong to only one shared group. Existing local linear/Darcy drag remains
active and converts kinetic energy to internal energy before thermal exchange.
The shared storage's own conductance is zero: every actual film belongs to its
declared interface.

Each group solves the simultaneous backward-Euler equations
`Q_i = dt G_i (T_wall,new - T_gas,i,new)` and
`Delta E_wall + sum(Q_i) = 0`. The unknown is an energy increment about the
initial wall state, not a difference of large absolute energy totals. Every
gas response reuses the original frozen-source caloric adapter and existing
safeguarded reservoir solver. Original source capacity integrals or T(h/u)
supply the wall storage and temperature, including latent plateaus; no fitted
gas Cv, apparent latent capacity, or repeated copy of the physical wall mass
is introduced. When original tangents exist, the scalar derivative uses the
individual gas heat responses and wall tangent; otherwise the same safeguarded
owner bisects. All interfaces retain the same solved endpoint wall state.
Reprojecting that state after convergence would change stiff interface fluxes
and is not performed.

Acceptance checks the actual represented wall and gas storage changes and
each original endpoint heat equation against the joule allowances. Large gas
formation-energy gauges cancel through the existing sensible-energy owner.
An original wall coordinate that cannot represent the required exchange must
meet the same actual-storage residual; a small solver update alone is never
acceptance. Original wall enthalpy retains its source's pressure-work convention
and is not relabelled as rigid-vessel internal energy.

The source uses the actual duration accepted by the transport retry.
Per-branch signed heat is added to the existing branch residual, while the
single wall inventory enters the global gas/wall and gravity-aware accounts
once. Fixed gas mass and volume produce no extra gravitational work during
heating. All branch/plenum fields, wall inventories, clocks and reports publish
together after the final cancellation poll. Source-domain refusal, insufficient
nested work, or cancellation preserves the preceding accepted history.
Prepared gas-response arrays and compensated reductions are reused.

`SharedWallControl` bounds the total original gas primitive/caloric and wall
state queries across all groups, brackets and final checks in one source
attempt; each original inversion retains its own bounds, and complete network
retries retain their existing cap. `SharedWallReport` retains the accepted
coupled storage residual, maximum individual wall and interface residuals,
actual source-query count, outer iterations and interface-solve count.

This is a fixed-geometry, spatially uniform common wall with explicitly
authored conductances. Transport/source splitting remains first order.
Initial wall temperatures must lie in the original gas caloric interval.
Spatial wall gradients, moving geometry, multiple shared groups on one gas
cell, inferred material conductivities and empirical heat-transfer qualification
remain separate capabilities. The scalar solve does not establish a general
discrete entropy theorem.

Seven focused physical regressions cover an independent three-body
backward-Euler solution with drag heating, nonlinear-Cv equilibration and stiff
passive response, large independent gas energy datums, exact hydrostatic and
thermal rest at a large common potential datum, one original latent inventory,
heat-driven field motion, source/work/overlap refusal, and final cancellation
with deterministic history replay.

### Shared walls with evolving mixture composition

With `experimental-euler`, the multicomponent network's
`wall_exchange::NetworkWallOperator::with_shared_walls` binds the same original
finite inventories and `SharedWallPort` values to evolving mixture cells and
plenums. The frozen and mixture adapters share one scalar energy-balance
kernel; every mixture query borrows the current constituent densities and
original gas caloric owner. Chemistry or transport between source halves
therefore changes the actual gas heat capacity and inversion used by the next
half, without freezing composition or replacing the source energy reference.

Local drag precedes shared heat. Every gas target must have an adiabatic local
contact without competing storage and may occur only once across the shared
groups. Any bound gas transport law also admits the heated endpoint. A moving
plenum uses the composing owner's current trial volume; only that owner can
change the accepted geometry, and its pressure work remains separate.

Candidate wall storage carries through both source halves and publishes with
the complete network transaction. Accepted accessors expose only the last
complete state. `SharedWallControl` counts original source queries across both
halves and resets on retry; `last_shared_wall_report` combines their storage
residuals, contact residuals and work. A failed second half cannot publish a
first-half wall state. Original h/u, latent energy, pressure-work conventions,
uniform-wall applicability, authored conductances and first-order composition
retain the limits described above.

## Original material conductivity on shared gas/wall ports

With `experimental-euler,material`, both the frozen and evolving-mixture
network wall owners expose `with_shared_wall_conduction`. Rows follow the
existing shared-wall order and entries follow each wall's gas-port order.
Each optional `MaterialWallConduction` retains its original temperature or
phase-dependent conductivity source and planar or radial resistance geometry.
A bound port's authored conductance remains the gas-side film; solid resistance is
added once in series. Unbound ports retain their complete authored coefficient.
The resulting conductance is never written back as the next step's gas film,
and no material density, heat capacity or latent energy becomes extra storage.

`InitialState` samples original k at the beginning of each source half and
checks its final state. `ImplicitEndpoint` resolves k(T(h),phi(h)) at each
proposed wall inventory inside the existing simultaneous gas/wall balance.
The phase form requires the same original equilibrium chart as the wall's
caloric owner, so conductivity can change on an isothermal latent plateau.
The existing temperature-coordinate form retains its original source and
admitted geometry. Implicit evaluation requires original caloric storage.

The shared energy-root kernel is unchanged in ownership and physical balance.
An implicit conductivity binding suppresses the incomplete capacity-only
Newton slope; the bounded scalar solver uses bisection without inventing dk/dh
or asserting nonlinear uniqueness. Every actual conductivity query spends the
existing source-attempt budget, including failed proposals and final checks.
If wall-to-wall links subsequently move the inventory, all bound conductivity
sources admit that final candidate before publication. Domain, work or
cancellation refusal preserves the previously accepted gas and wall states.

This is a uniform-wall resistance model, not a resolved through-wall thermal
profile. The material supplies conductivity and the caller supplies the actual
resistance geometry and gas film. Source support at sampled states does not
establish a continuous-path certificate or empirical heat-transfer accuracy.

## Connected finite wall inventories

`NetworkWallOperator::with_shared_wall_links` connects the existing shared-wall
inventories through ordered `SharedWallLink` thermal ports. Each has two
distinct original wall indices and one nonnegative authored conductance in
W/K, complete unless material resistance is explicitly bound below. Cycles and
parallel interfaces retain authored order; a zero-conductance
link transports no heat. No gas mass or species passes through these ports.

After the original simultaneous gas/shared-wall solves, each link solves
`Q = dt G (Ta(-Q) - Tb(Q))` with the existing safeguarded scalar owner. Original
capacity integrals or T(h/u) charts supply both endpoint states, including
latent plateaus and their conserved phase energy. Exact available tangents
accelerate the same equation; absent local tangents use its bounded bisection.
For fixed conductance the heat interval intersects both original source domains
and the passive initial-flux bound. No source extrapolation, fitted capacity, duplicate wall
mass or alternate nonlinear solver is introduced.

Acceptance checks both actual represented storage changes and the stage heat
equation. Their signed sum enters the existing whole-network wall-energy
account; wall-to-wall transfer is internal and never counted as external heat
or gas heating. Source queries share `SharedWallControl` across gas groups and
all link stages, including failed evaluations and final checks. Prepared
receipt arrays are reused. Failure or final cancellation leaves gas fields,
wall inventories, clocks and all previously accepted reports unchanged.

`shared_wall_link_reports` retains every link's initial and endpoint original
states and temperatures, signed heat received by the second wall, both actual
storage changes, energy/flux residuals and scalar iteration count. These are
stage endpoints: a later link can change a wall again. The final wall state is
not substituted into an earlier gas or wall heat equation. All links use the
actual interval accepted by the original network transport/retry owner.

This is first-order composition of finite thermal ports, with observable
link-order error that must be controlled by timestep refinement. It is not a
simultaneous backward-Euler solve of the whole thermal graph, nor a spatial
PDE discretization. Geometric/material spatial conduction remains owned by
`fs-thermal`; no conductivity or geometric resolution is inferred from the
authored W/K ports. Existing fixed-geometry and lumped-cell validity conditions
remain. The scalar solves alone do not establish a general entropy theorem.

G0/G1 regressions cover independent staged gas/wall calorimetry; an analytical
three-wall chain with first-order convergence and vanishing link-order bias;
latent inventory transfer and exit into sensible cooling; exact thermal rest
with cyclic/parallel links and isolated storage; exact whole-attempt source
budget boundaries and late source-domain refusal; and atomic cancellation with
bitwise checkpoint replay. Source-specific calorimetry, energy closure and
temperature-hull checks do not claim experimental heat-transfer accuracy.

### Connected walls in evolving-mixture transactions

The multicomponent network wall owner also exposes `with_shared_wall_links`,
`shared_wall_links` and `shared_wall_link_reports`, using the same port types
and scalar exchange kernel as the frozen network. Each source half first
solves its current gas/shared-wall groups, then the authored wall links. The
next half starts from those candidate wall states and the newly transported
or reacted gas composition. Wall-to-wall heat changes only the two original
wall stores; it adds neither external heat nor another gas heating term.

The complete source-attempt budget covers both gas and wall links across both
halves. Report buffers are allocated at binding for two full link sequences.
The accepted report slice contains one chunk in authored link order for each
source half: one chunk for standalone `advance_for`, two for a composed
network step. Every receipt retains its actual half duration and original
stage endpoints. A later stage cannot replace the temperatures of an earlier
flux equation. Accepted states and the full report sequence publish together;
failure or cancellation preserves the previous complete transaction, and a
retry resets every candidate state, receipt and source-work counter.

Focused tests check independent three-wall staged calorimetry, transfer out of
an original latent plateau into sensible cooling, both-half receipts and exact
source-work limits, final cancellation with replay, and source-domain refusal
after gas heating. The same first-order, uniform-wall, authored-conductance and
source pressure-work limits apply. This extension adds no spatial gradient,
empirical heat-transfer qualification or simultaneous graph-endpoint claim.

### Original material resistance between finite walls

Both network wall facades expose `with_shared_wall_link_conduction` and
`shared_wall_link_conduction`. Rows follow authored link order; each row's
two optional `MaterialWallConduction` bindings belong to that link's first
and second wall. Each side retains its own original caloric chart, phase
state, conductivity source and length/area. Distinct materials never share an
averaged temperature, enthalpy datum or inferred conductivity law.

On a bound link, the original authored G represents finite interface/contact
conductance. Each supplied bulk resistance L/(k A) is added once in series
with it using scaled conductance arithmetic. The authored G remains unchanged;
zero still means insulation, and two absent layers retain the complete-G path.
The finite interface model does not infer perfect contact. Geometry describes
each bulk segment between its retained cell state and the common interface;
it creates no additional mass or latent-energy account.

Initial-state layers sample at that link's actual entering state, including
changes made by earlier links. Implicit layers instead query each proposed
endpoint; the two sides may use different evaluation modes. When conductivity
is implicit, the heat bracket uses original passive storage limits rather than
the initial-flux cap, because increasing k can make the accepted heat exceed
that initial estimate. Safeguarded bisection avoids inventing dk/dh or claiming
nonlinear uniqueness. Fixed and lagged links retain their available tangents.

Every runtime conductivity query shares the complete gas/wall source budget.
Lagged endpoints are admitted, and all bound link sources are checked again
after the full ordered sequence: a later link cannot leave an earlier source's
wall outside its support. Both facades retain the same original storage,
physical heat receipts and atomic publication. The first-order composition,
fixed-geometry, finite-cell and material-domain limits above still apply.

Focused G1/G3 cases compare independent two-material calorimetry in both
orientations and an increasing-conductivity latent root beyond the initial-flux
bound. G4 cases check source-work limits, cancellation/replay and late-link
material-domain refusal in both facades. Actual-source native compilation and
28 focused gas tests pass; workspace/DSR verification remains separate.

## Finite-time thermal exchange with multicomponent phase equilibrium

With `material,experimental-euler`, `material::phase_thermal::exchange_mixture_parcel`
couples the original `IsochoricMixtureEquilibriumParcel` to a held-temperature
reservoir, an authored constant-capacity finite wall, or one retained
`MaterialWallCaloric`. The parcel's immutable feed, mass and vessel volume
remain unchanged; pressure, temperature and phase amounts follow the original
source-bound U,V equilibrium. Every trial uses that same caloric owner and its
finite phase-screen evidence. No equilibrium Cp/Cv, dilute transport, extra
latent heat, chemical reaction or phase-change kinetics is inferred.

`exchange_isobaric_mixture_parcel` uses this same implicit thermal owner with
the original `IsobaricMixtureEquilibriumParcel` and `IsobaricMixtureThermalControl`.
Its PH source resolves heat as Delta H, including latent phase redistribution,
while volume changes against the prescribed constant pressure load. The source
independently supplies Delta U and actual outward `p Delta V`; the final joint
gas/wall first law includes that work exactly once. The common heat report's
mechanical work to gas is `-p Delta V`, while `heat_to_gas_j` is Delta H.
The original isobaric caloric receipt exposes both energy changes and work.
Finite-wall exchange remains equal and opposite even through a latent plateau.
Mass, species and pressure stay fixed; changing contact geometry, piston inertia
and evaporation kinetics remain outside this quasistatic lumped model.

The existing bounded wall scalar owner solves
`Q = dt G (Twall(-Q) - Tgas(Q))` in signed heat coordinates. `G` is the complete
authored conductance in W/K, fixed during the supplied physical interval.
The caller declares a heat bracket; it is intersected with the initial passive
flux bound and original finite-wall storage limits. Every remaining gas endpoint
must satisfy its original pressure/temperature/source brackets. An unsupported
trial, missing sign bracket or unresolved original phase search refuses rather
than extrapolating or switching to a homogeneous approximation. Numerical
secants accelerate the safeguarded scalar root; they carry no heat-capacity or
physical-tangent authority. The accepted source-bearing trial is retained and
published without another caloric solve.

Actual represented gas energy, original finite-wall h/u changes and endpoint
heat flux must independently close in joules. Original wall enthalpy keeps its
declared pressure-work convention; it is not relabelled rigid-volume internal
energy. Reservoir heat is external, finite-wall exchange is internal, and
latent storage already belongs to the original thermodynamic potentials.
The U,V or PH residual allowance reserves part of the thermal budget; a nonzero
initial flux tightens the absolute root allowance so a small exchange does not
vanish into a loose configured tolerance. Source increments cancel species
formation and entropy datums before rounding.

`MixtureThermalControl` bounds heat queries and all original PT selections for
the complete interval, in addition to the existing per-trial pressure/phase
limits. Remaining PT work is passed into every original caloric solve. Reports
retain actual aggregate heat/PT/pressure work and the sum of original EOS work
ceilings, separately from the accepted trial's caloric report. Final
cancellation, source/balance failures and work exhaustion leave both supplied
inventories unchanged. Zero duration, zero conductance and equal initial
temperatures retain both states exactly with zero phase evaluations. A cloned
parcel and retained original wall state supply deterministic retry state.

This is a first-order backward-Euler, spatially uniform thermal interface.
Uniform gas/wall assumptions, authored conductance applicability and time
refinement remain caller obligations. It does not establish spatial conduction,
evaporation rate, slip, dynamically resolved boundaries, global Gibbs stability, calibrated
heat-transfer correlations or experimental fluid accuracy. The gas source
entropy increment is retained; an arbitrary finite wall may lack its own
entropy law, so the API makes no combined entropy-production theorem.

The public `phase_equilibrium_thermal` integration target adds six focused
G0/G1/G3/G4 tests: independent binary isochoric calorics with a held bath;
finite-wall heating across the phase boundary; original latent-wall cooling;
independent species-datum and extensive-scaling invariance; exact no-work rest
and finite source-domain refusal; and whole-call work, final cancellation and
deterministic retry. The methane/ethane fixture uses synthetic NASA capacities;
its original cubic/fugacity references test constitutive consistency, not a
measured-mixture validation claim.

Four additional isobaric regression groups verify independent bubble/dew
transition references, reservoir and finite-wall heat, Delta H/Delta U/actual
pressure work, reversal into a latent wall, large species-datum invariance,
source/work refusal, exact rest and final-cancellation replay. They exercise the
public consumer and shared scalar owner, with no second phase or heat model.

Native validation compiled the production path and passed five of these six
cases. The finite-wall phase-crossing fixture now declares a supported heat
bracket and a 256-iteration original flash budget, retaining the original
strict tolerances. Its final rerun was blocked by shared disk exhaustion and
remains unverified; neither that case nor the complete suite is claimed green.

## Heat exchange between finite phase-equilibrium mixtures

`material::phase_thermal` also couples two finite original mixture inventories
through `exchange_isochoric_mixture_pair`, `exchange_isobaric_mixture_pair`, and
`exchange_isochoric_isobaric_mixture_pair`. One implicit heat coordinate closes
`Q = dt G (T_right(-Q) - T_left(Q))`. Positive Q heats the left participant;
each side keeps its own species, pressure or volume constraint, original phase
policy, and source calorics. Heat changes U for a rigid participant and H for a
constant-pressure participant. Latent storage stays inside those original laws.

`MixturePairThermalControl` declares independent source bounds and one shared
heat/PT work budget. Both actual source heat increments, both endpoint heat
equations, combined Delta U plus both outward pressure-volume works, and the
sum of original entropy increments are checked before either parcel changes.
The entropy tolerance permits only a declared numerical deficit. Neither
external heat nor mass/species exchange enters this isolated thermal pair.
Each prescribed pressure load receives its own boundary work.

The last admitted pair is retained without an additional equilibrium solve.
Zero duration/conductance or equal temperatures preserve both states without
PT work. Source refusal, budget exhaustion, balance refusal, and final
cancellation preserve both original inventories. This backward-Euler interface
assumes uniform inventories and held conductance, with time refinement required;
it does not supply spatial heat-exchanger geometry or phase-change kinetics.

## Steady heat exchange between perfectly mixed mixture streams

With `material`, `material::phase_heat_exchanger::exchange_steady_mixture_heat`
joins two steady, perfectly mixed chambers through an authored total thermal
conductance `UA` [W/K]. Each chamber retains its own original mixture source,
fixed pressure, positive mass flow and phase policy. Positive duty `Q` enters
the left stream. The existing pressure/enthalpy increment owners resolve
`Delta h_left = Q / mass_flow_left` and
`Delta h_right = -Q / mass_flow_right`; a bounded safeguarded scalar solve
closes `Q = UA (T_right_out - T_left_out)`.

The original NASA increments, EOS departures and changing phase amounts carry
sensible and latent energy together. Independent species enthalpy and entropy
datums cancel before rounding. Stream enthalpy includes pressure-flow work;
there is no additional `p Delta v` debit, fitted mixed heat capacity, or latent
heat source. The returned original PH receipts preserve both source laws and
the finite phase-screen evidence. Actual `mass_flow * Delta h` from each source
must close both interface equations and the combined energy balance. The sum
of original `mass_flow * Delta s` must be nonnegative within a declared numerical
allowance and is reported without clipping. This total includes irreversibility
from mixing each inlet into its chamber as well as interface heat transfer.

The signed duty bracket must contain zero. The passive initial-flux bound
further restricts it, and both resulting endpoints must lie within the original
source domains. Native PH tolerances are tightened against the power allowance
and each stream's mass flow. A shared duty-query limit and total PT allowance
cover both streams, the once-prepared inlets and all bracket/root trials; each
query receives only the remaining work allowance. An immutable prepared PH path
reuses each inlet's original phase selection without changing its source or
phase policy; its one-time work is included in the exchanger totals. Native source/domain refusal,
failed balances, exhausted work or final cancellation returns no outlet pair.
The accepted trial is retained without another flash. Zero conductance or equal
inlet temperatures returns exact zero duty after one original PT selection per
stream. Mass flow remains strictly positive even at rest.

This model assumes no shaft work, pressure drop, external heat, wall storage,
or change in stream kinetic/potential energy. It describes two mixed chambers;
it does not describe counterflow or parallel-flow plug-flow geometry, predict
conductance or flow, introduce phase kinetics, or certify global phase stability.
Focused regressions use independently calculated binary PR phase/caloric
references, a liquid-to-coexistence transition, the analytic dilute constant-Cp
limit, unequal and reversed streams, independent large reference-datum shifts,
exact rest, shared work boundaries, domain refusal and cancellation/replay.

## Spatial parallel and counterflow mixture heat exchangers

With `material`, `material::phase_heat_exchanger::spatial` resolves the
one-dimensional steady plug-flow energy balance for parallel or counterflow
streams. Each stream retains its own mixture, fixed pressure, positive mass
flow and original phase policy. Two immutable prepared PH paths select the
inlets once. All local states use native source-integrated increments from
those same inlets, including NASA calorics, EOS departures and phase conversion.

For signed total heat duty `Q` into the left stream and cumulative heat fraction
`u`, the left increment is `u Q / mdot_left`. The right increment is
`-u Q / mdot_right` in parallel flow and `-(1-u) Q / mdot_right` in counterflow.
The continuum relation is `UA = Q integral[du / (T_right(u) - T_left(u))]`.
Thus each species feed remains separate and both stream energy balances close
without chaining approximate outlet states, fitting a heat capacity, adding a
latent source, or counting pressure-flow work twice. Counterflow outlet
temperatures may cross; passivity concerns the co-located temperature gap.

The bounded scalar root uses the power residual `Q - UA / I`, with sign
normalization and `I` the positive reciprocal-gap integral. Adaptive nested
Simpson panels cache original endpoint and midpoint states. Native PH work is
linear in sampled nodes per duty query. Breadth-first subdivision avoids
exhausting one near-pinch interval before the complete partition can exclude an
excessive trial duty. Monotone source-temperature bounds require each accepted
panel's entire gap to be positive under the native monotone caloric inversion
assumption. A partition-wide lower bound on required conductance can reject an
upper duty without fully resolving its near-pinch profile. Witnessed
nonpositive gaps also reduce only the upper bracket. Original source/domain
failures and unresolved panel bounds remain explicit refusals.

The raw coarse/fine quadrature difference is an estimate, without a smooth-order
Richardson divisor that could hide phase-boundary kinks. Caller conductance
tolerances are further limited by the reserved power budget. Acceptance requires
the scalar residual plus its conditional quadrature power uncertainty to fit
the physical allowance. This is not a certified quadrature bound and excludes
native PH inversion uncertainty and model/source-data error. Original PH
tolerances remain explicit. Actual source enthalpy flows and combined source
entropy are checked by the same final balance owner as the mixed-chamber API.

The returned source-bearing profile distinguishes cumulative heat fraction from
cumulative conductance fraction. The latter equals axial fraction only for
uniform conductance per unit length. Its normalization uses the same
compensated reduction as the accepted total integral. Both physical ends,
original phase states, source increments, total PT/EOS work, quadrature estimates
and minimum sampled gap are retained without another accepted-state flash.

A single whole-call PT budget includes both preparations and every duty/root
sample. Separate duty, panel and refinement-depth limits bound the remaining
work and storage. Strict monotone gap bounds can require increasing panel counts
at high NTU; an insufficient declared budget refuses rather than silently
skipping pinch checks. Zero UA or equal inlet temperatures gives exact zero
duty after the two inlet preparations, with no new endpoint PT calls.
Cancellation or any numerical/source refusal publishes no partial profile.

The model includes no pressure drop, axial conduction, wall heat capacity,
external heat, shaft work, species exchange, reaction kinetics or change in
kinetic/potential energy. It retains the original finite phase-screen evidence,
without a global phase-stability claim. Focused tests cover independently
calculated PR parallel/counterflow references, a liquid-to-coexistence transition,
the dilute-gas effectiveness-NTU formulas, high-NTU counterflow outlet crossing
and near-pinch upper bounds, source-reference shifts, reversed flow duty, exact
rest, shared work limits, spatial refinement, source refusal and cancellation.

## Adiabatic equilibrium mixture compressors and expanders

With `material`, `material::phase_machinery::adiabatic_mixture_machine`
resolves a steady adiabatic compressor or expander over one immutable
`ThermochemicalPengRobinsonModel`. It first solves the original pressure/entropy
increment constraint `Delta s = 0` at the outlet pressure. For an authored
isentropic efficiency `0 < eta <= 1`, the actual enthalpy target is
`Delta h_s / eta` for compression and `eta * Delta h_s` for expansion. The
existing pressure/enthalpy increment owner resolves the actual outlet using
exactly the same original phase-screen policy. Source intervals may differ
between the ideal and actual temperature solves; phase policies may not.

Both solves retain the complete feed, original NASA laws, EOS departures,
phase mixing entropy and finite original phase-screen evidence. Liquid/vapor
amounts may change. Species energy and entropy reference datums cancel through
the original integrated increments before rounding. No mixed heat capacity,
extra latent energy, dilute transport extrapolation or replacement EOS enters
the machine model.

The signed shaft work **into** the fluid is the actual source-integrated
enthalpy rise. This steady-flow balance assumes no external heat and negligible
inlet/outlet kinetic and potential energy differences. Stream enthalpy already
contains pressure-flow work; another pressure-volume debit would double count
it. The result retains the ideal endpoint, actual endpoint, source work,
actual entropy production, enthalpy residual and ideal entropy residual.
The caller may convert work to power with an independently supplied nonnegative
typed mass flow. Pressure ratio and efficiency do not predict that flow.

Admission checks the pressure direction and representable work sign. Actual
source entropy change must be nonnegative within a declared numerical
allowance; the reported value is never clipped. A single bounded PT allowance
covers both original solves, including repeated inlet queries and all bracket
states. Each nested owner receives only the remaining allowance. Cancellation,
source-domain refusal, unresolved phase work or an exhausted budget returns no
partial outlet. Efficiency one reuses the accepted PS endpoint; equal pressure
returns exact zero work and entropy change after one original inlet selection
when the reference is inside its admitted bracket.

This is an equilibrium endpoint model with a caller-authored efficiency, not a
validated machine map, rotor/transient model, reaction/phase-kinetic model,
resolved flow geometry, choking model, or global phase-stability certificate.
The declared tolerances bound numerical residuals, not model-form error.

Seven focused regressions cover an independent two-phase isentropic expansion,
reversible round trip and both efficiency conventions, actual nonnegative
entropy production, large independent species reference shifts, exact rest,
whole-call budget boundaries, malformed policies and final cancellation/replay.
Rustfmt parsed the new source and tests. Native compilation and execution remain
unverified because the shared build disk is exhausted; no passing native result
is claimed for this new module.

## Source-bound real-gas normal-shock compression

`material::nonideal_shock::normal_shock` joins the original Peng--Robinson
mechanics, NASA calorics and material contexts to the normal-shock mass,
momentum and energy equations. It solves the Rayleigh pressure against the
native density/internal-energy reconstruction on an explicit compression
bracket. The Hugoniot energy uses the original finite internal-energy increment,
including both density departures; independent species formation datums cancel
before rounding. Native sound speed, Cv and thermal pressure derivatives supply
the safeguarded Newton tangent and sonic checks. No fitted gamma enters.

Both density bounds exceed the inlet density, and a small momentum residual
cannot admit the identity root: upstream flow must be supersonic and downstream
flow subsonic in the shock-relative frame. The complete candidate must pass
mass, momentum, Hugoniot energy, stagnation enthalpy and signed entropy gates.
Numerical tolerances, temperature/density bounds and iteration work are explicit.
Source refusals or cancellation publish no candidate and preserve all inputs.

The original composition, critical properties, pair interactions and material
context remain attached to both resolved states. Five focused regressions cover
an independent high-precision dense PR shock, the analytic dilute limit,
formation/entropy reference invariance, source/bracket/work refusals and
cancellation/retry. This assumes one caller-declared continuous homogeneous
stable branch; local stability checks do not establish global phase stability.
It does not resolve a two-phase front, chemical or vibrational relaxation,
shock thickness or spatial field evolution.

### Real-gas compression from prescribed piston motion

`material::nonideal_shock::piston_compression` instead takes signed upstream
and piston velocities along a normal directed from the piston into the gas.
The piston must move toward the upstream gas in that relative coordinate.
The original bounded density solve now determines the shock speed and wall
pressure together; the caller does not supply a guessed shock Mach number.
Each trial uses the exact piston mass/momentum relation and finite Hugoniot
energy increment, with an analytic tangent from the same native PR calorics.
The complete shock retains the original source states, conservation defects
and entropy/sonic admission, plus the laboratory velocities and shock speed.
Its embedded normal-shock speeds and mass flux are positive magnitudes;
the moving-shock mass flux along the declared laboratory normal is negative.
The known velocity increment is retained in momentum and kinetic differences
to avoid subtracting nearly equal large shock-relative speeds.

This is the compressive piston branch on the same declared homogeneous
source domain. It does not supply an expansion fan, advance a moving mesh,
integrate piston inertia or evolve a spatial real-gas field. A common change
of laboratory velocity leaves the physical jump unchanged within represented
input precision; unrepresentable shock/wall separation is refused.
Additional authored regressions cover a separate 75-digit dense PR piston oracle,
the analytic dilute piston limit, changes of reference frame, and bounded
source/branch/work/cancellation refusal. These numerical reference fixtures
do not establish empirical calibration of the authored material mixture.

## Simultaneous reaction and finite reservoir heat exchange

The prepared constant-volume and constant-pressure reactors accept an explicit
`HeatBath`: nonnegative total conductance G [W/K] and a held positive bath
temperature Tb. The original extent integrator advances Qdot=G(Tb-T) together
with the reaction laws. Both backward Euler and SDIRK2 eliminate stage heat
inside the original caloric root; constrained temperature derivatives use
original reaction delta-u/(Cv+hG) or delta-h/(Cp+hG), with reference-volume
scaling. Physical heat-capacity reports remain Cv/Cp. There is no additional
integrator, split temperature update, reaction heat source or heap workspace.

Original initial U/H targets retain their meanings. Accepted energy balances
subtract compensated external heat once; fixed-pressure results retain actual
volume and p delta-V work. Full, half and rejected stages carry independent
heat state. Step doubling checks heat as well as species and temperature using
a formation-datum-independent scale. Caloric residuals and local heat
differences retain the compensation terms.

Parcel entropy may decrease on cooling. Admission uses parcel-plus-reservoir
production delta-S-Q/Tb; reports retain parcel delta-S, reservoir transfer,
total production and minimum accepted-step production. Original element,
charge, mass, nonnegative inventory, source support and energy gates remain.
Zero conductance with radiation disabled preserves adiabatic results.
Explicit-energy zero-duration calls resolve initial U/H before any exchange.

Six new native tests passed: frozen heating/cooling, simultaneous reaction and
heat against an independent closed form including resonance, fixed-pressure
first law/work, large elemental formation-reference shifts, insulated parity,
cancellation/replay, explicit targets and constrained-Jacobian differences.
The bounded native fs-gas library run passed 68 tests, with two lengthy existing
hydrogen/methane reference cases excluded; the SDIRK hydrogen reference passed.
These are focused G0/G1/G3/G4 checks, not a full repository Gauntlet result.

`HeatBath::with_gray_radiation(area, effective_emissivity)` additionally admits
gray thermal radiation to the same held reservoir. The caller supplies a
nonnegative area and an effective emissivity in [0,1], including the declared
gray geometry closure. The bath retains their product epsilon*A [m²]; it does
not infer opacity from the reacting species. The complete exchange is
Qdot=G(Tb-T)+sigma*epsilon*A*(Tb^4-T^4), using the same rounded SI constant
sigma=5.670374419e-8 W m^-2 K^-4 as the native gray surface/gap owners.

Both original implicit methods eliminate that complete stage heat inside the
same total-U or total-H root. Its analytic caloric feedback is
dt*(G+4*sigma*epsilon*A*T^3), with the original reference-volume scaling.
Fourth-power differences are factored by the represented Kelvin separation;
no large-power subtraction, separate radiation timestep, frozen radiation
conductance, or duplicate formation-energy source is introduced. Radiation
and conduction share the original compensated heat and held-reservoir entropy
account. Zero effective area retains the original linear arithmetic exactly.
The declared area remains fixed during an advance even if the gas expands.
Original material parcel and piston adapters forward the same complete bath.
`HeatBathSchedule` linearly interpolates the physical epsilon*A coefficient
alongside G and Tb, retaining authored baths at exact knots; it does not infer
independent area and emissivity trajectories.

Added native regression cases cover the independent analytic fourth-power
cooling/heating integral, simultaneous reaction/radiation against a separate
scalar ODE reference, both U/H constraints and time methods, analytic
constrained-Jacobian differences, near-equilibrium fluxes, admission, and a
programmed radiating-area material parcel with heat/work/entropy closure.
These numerical fixtures do not establish empirical radiation calibration.

Scope remains a sealed homogeneous ideal-gas parcel and an infinite
fixed-temperature reservoir, with caller-declared conductance and optional
lumped gray exchange. This bath does not imply finite wall capacity,
participating-gas spectral transport, an evolving wall source, streams,
spatial transport or empirical mechanism validation. Time-error estimates are
local. Finite original-source bodies remain the separate owner below.

## Finite original-mixture gas thermal partners

`chemistry::heat::gas::FiniteGasHeatBody` borrows one original
`FrozenIdealGasMixtureModelV1`, fixed molar amount, initial temperature and
independently admitted temperature interval. `GasHeatConstraint` chooses rigid
volume or held pressure for this partner, independently of the reacting gas.
Both prepared reactors expose `advance_with_finite_gas_body`; the original
reactor and continued partner endpoints publish together. Partner composition
is frozen, no species cross the interface, and declared conductance/optional
gray epsilon*A have the same explicit geometry convention as the held bath.

The original NASA-mixture delta-u/delta-h and Cp/Cv enter the existing
backward-Euler/SDIRK2 stages. Source formation references cancel before
arithmetic. One bounded, cancellable monotone caloric root eliminates partner
temperature with the analytic nonlinear heat feedback; the original adaptive
integrator retains cumulative heat, trial rejection and total entropy admission.
Partner-temperature error is included in full/half comparisons. Caloric roots
must resolve both the existing Joule tolerance and a share of the existing
relative-temperature allowance, so a small finite inventory cannot disappear
inside a larger reactor's energy budget. No per-trial allocation is introduced.

Rigid partner storage satisfies delta-U=-Qgas. Held-pressure storage satisfies
delta-H=-Qgas and separately retains delta-U and outward pressure work. Work
uses n*R*delta-T, avoiding subtraction of rounded endpoint pV. Original
Cv/T or Cp/T entropy differences retain the partner's actual temperature path.
Reports include heat, H/U changes, work, entropy and heat residual; original
amount, mass and source identity remain unchanged. Source positivity and
caloric continuity are admitted across the complete partner interval. The gas
temperature may lie outside that interval while the actual partner remains
inside it; source exhaustion refuses before publishing either endpoint.

Focused regressions cover all four mechanical-constraint combinations,
reacting/frozen two-body closed forms, nonlinear source Cp and radiation,
analytic feedback, large energy/entropy reference shifts, independent source
domains, tiny-inventory accuracy, cancellation/replay and zero-coupling parity.
This is a homogeneous ideal-gas thermal partner, with frozen composition and
fixed declared coupling geometry. It implies no reacting second mechanism,
mass exchange, spectral opacity or empirical heat-transfer validation.

## Finite material bodies in adaptive reacting parcels

With the `material` feature, `chemistry::heat::material::FiniteHeatBody`
retains an exact original equilibrium enthalpy chart, invariant body mass,
declared constant body pressure and total thermal conductance. Both prepared
gas reactors expose `advance_with_finite_body`; the returned gas and original-
source body endpoints publish together. Source, cancellation, work and energy
refusals publish neither endpoint.

Without radiation, the existing backward-Euler/SDIRK2 chemical stages eliminate
the finite body on its original caloric intervals. Retained affine Cp(T) uses a
scaled quadratic solve of h(T)+beta*T; authored affine T(h) and latent intervals
retain their linear elimination. This original analytic path is unchanged
when radiation is disabled. The constrained
gas Jacobian uses the original derivative at the accepted body temperature,
including source Cp variation within compact two-knot charts. Latent plateaus
retain their actual enthalpy and phase fraction; the adaptive comparison
includes both body temperature and liquid
fraction. Caloric bracketing respects both original source domains.

`FiniteHeatBody::with_gray_radiation(area, effective_emissivity)` adds the same
explicit gray geometry/emissivity convention as the held bath, retaining
epsilon*A alongside the original body state. The evolving body exchanges
G*(Tb-Tg)+sigma*epsilon*A*(Tb^4-Tg^4) with the reacting gas. Its temperature
responds to the transferred heat, including finite sensible storage and sharp
latent storage. Fixed radiation geometry survives continued body states and
prescribed conductance histories; density changes do not silently rescale area.

The same chemical stage retains cumulative gas heat q and body
h=h0-V*q/m. Binary lookup locates the original caloric segment. A bounded,
cancellable safeguarded Newton solve uses heat relative to the fixed stage
history, evaluating the original h(q) source and complete gray flux together.
Sharp latent segments retain direct elimination. The exact eliminated gas
feedback is `(a+4*b*Tg^3)/(1+(V/m)*(a+4*b*Tb^3)*dTb/dh)`, with
`a=dt*G/V` and `b=dt*sigma*epsilon*A/V`. Original affine Cp(T) and one-sided
source derivatives remain authoritative; no apparent latent heat capacity is
introduced. Monotone endpoint solves constrain the outer gas-temperature
bracket to the original body domain. All new roots use the caller's existing
energy tolerance and maximum caloric iterations, and perform no allocation.
This adds no separate heat timestep or duplicated radiation-energy ledger.

Gas heat and body heat cancel once. The original material isobaric heat owner
reports density-derived volume, pressure work, enthalpy and internal-energy
changes. Fixed-volume gas closes its U change with body U and body p*dV;
fixed-pressure gas also reports its own outward p*dV. Formation and latent
energy remain in their respective original potentials. Entropy admission uses
the actual gas change plus the original body integral of dh/T. Unresolved
enthalpy-datum heat increments refuse rather than disappear.

This is a homogeneous finite-body model with declared conductance, optional
gray exchange, held pressures and no gas/body mass exchange. The caller must justify negligible
internal body gradients and the pressure applicability of the original chart.
It does not resolve wall geometry, contact mechanics, moving interfaces or
spatial flame transport. Source traversal is bounded to 4096 original knots.

Five focused G1/G3/G4 tests cover independent sensible/reacting two-body
solutions, sharp latent heating/cooling and onset, both gas constraints and
time methods, enthalpy datum changes, both pressure-work accounts, finite-body
constrained Jacobian differences, source refusal, cancellation, replay and
insulated agreement with the original adiabatic owner. The final actual-source
material-feature heat target passed all 11 tests (five new, six existing).
A broader material-feature unit run passed 78 tests, excluding the two known
long hydrogen/methane Cantera cases, before the final finite-denominator guard
and additional Jacobian test; that final focused run includes both changes.
Default-feature compilation and all six existing heat tests also passed.
These are native fallback builds against compatible cached dependencies;
DSR and required sibling checkouts were unavailable. They do not establish
a full workspace verdict or measured-material validation.

The `source_cp_body` regressions additionally cover increasing/decreasing Cp,
original-interval heat closure and finite-difference feedback, shifted enthalpy
datums, and independent two-body evolution through both reactor constraints and
time methods. The focused original-source material-feature heat build now
passes all 21 tests, including these Cp regressions and the four finite-gray
cases. A separate default-feature heat build passes all 10 tests. The material
build retains original phase caloric/entropy and chemistry code while omitting
unrelated adapters. Full material/workspace Cargo execution remains blocked by
the missing sibling asupersync checkout.

Additional finite-gray regressions compare source-Cp stage heat and analytic
feedback against independent polynomial/finite-difference calculations, and
compare reacting/frozen two-body evolution against independent scalar ODEs.
They include sensible-to-latent heating and freezing, both gas constraints and
time methods, pressure work, original body entropy, source endpoint brackets,
zero-radiation parity, explicit root work limits, cancellation and replay.
These numerical fixtures do not validate a measured surface emissivity or
justify the homogeneous-body approximation for an arbitrary specimen.

### Finite support inertia in the original network wall transaction

`NetworkWallOperator::with_finite_supports` binds optional independent
`SlidingSupport` masses/velocities in the original branch/cell order. Selected
supports replace prescribed wall motion, retain authored or source-refreshed
resistance, and use the original exact reduced-mass recoil mechanics. By default
all slip loss thermalizes once in gas sensible energy. The existing local,
material-caloric, shared-thermal and conduction owner then advances heat with
drag disabled; no thermal inventory is copied or reconstructed. Both material
transport endpoint checks use actual support-relative slip, and the accepted
contact retains the returned support velocity.

`with_finite_support_heat_partition` instead supplies each independent support's
declared gas fraction. The remainder heats its original finite local wall before
the existing conduction/heat owner runs. Material capacity uses its exact
integral and inverse; native internal-energy plateaus retain latent inventory.
Mechanical closure admits internal energy or explicitly negligible expansion,
not unqualified enthalpy or generic sensible capacity. The wall's thermal amount
is already part of the support inertia. Shared thermal inventories require their
own coupled owner and cannot be duplicated as local frictional heat stores.
The thermal phase retains the original film and conduction owner, including
implicit endpoint evaluation when selected. Its receipt uses the complete original-to-returned wall
storage difference; frictional wall heat is an internal transfer out of gas,
not an additional source. All-gas partition preserves the existing default.

Support candidates follow the same begin-trial, source-half, retry and publish
lifecycle as gas and walls. Workspaces allocate at attachment; stepping adds no
allocations or clock. Read-only support states and per-cell `SupportExchange`
receipts retain actual complete support endpoint impulse/kinetic change and
accumulated friction-only gas impulse/energy and wall heating. Its friction-stage
energy balance includes gas, support kinetic and wall storage. Independent
dimensional momentum and energy budgets reject before publication.
`NetworkWallReport.mechanical_work_to_gas_j` includes finite-support work as
internal exchange with that support; the support receipt closes this additional
kinetic storage. Branch axes are not summed.

Each selected cell has its own mechanical support; identical descriptors do not
share inertia. Shared thermal storage remains separate. The feature does not
move geometry, evolve plenum bulk momentum, infer friction coefficients or
supply multi-cell rigid support dynamics. Seven existing selected-source native
tests cover analytic recoil/heat balance, both halves, cancellation/retry,
material-relative domain checks/Galilean boosts, identity, and shared thermal
composition. Four additional tests cover nonlinear wall heat partitions against
an independent quadratic heat balance, latent storage, prior accepted-state
rollback/replay, and physical admission. All eleven passed with the original
network and mixture implementations; optional material-feature and full Cargo
verification remain separate.

### One finite support exchanging momentum with multiple gas amounts

`SharedRecoilWorkspace` advances fixed gas amounts against one common
`SlidingSupport` along an authored scalar axis. Each `SharedRecoilCell`
retains its own density, volume, momentum and nonnegative frozen linear and
quadratic resistance coefficients. For unequal drag laws or heterogeneous nonlinear states, the original exact
reduced-mass pair mechanics is composed in a palindromic Strang sweep. This general path
approximates simultaneous shared-support dynamics with one common inertia.
Single-active-cell intervals reuse the exact original pair map. Multiple
active cells with one common linear rate and no quadratic drag use the exact
collective/deviation modes described below.

Bounded step doubling estimates local errors in each velocity, physical slip
heat and Galilean-invariant work. The support velocity is included. Dimensional
velocity/energy tolerances and a relative tolerance are separate from the
strict endpoint momentum and energy budgets in `RecoilControl`. On the general
refined path, each candidate interval also resolves the resistance scale:
`h*(max(k_i) + sum(m_i/M*k_i)) <= 0.5`, with
`k_i=a_i+2*b_i*R` and R the initial velocity-hull width of the support and active cells. This prevents
saturated split pair maps from masquerading as resolved heat allocation.
Rejected attempts consume the explicit step budget; under-resolved,
non-advancing or exhausted intervals refuse. These are local time-error
estimates, not a certified global or stiffness-independent accuracy bound.

Every exact pair flow exchanges opposite impulse and dissipates nonnegative
physical slip energy. Accepted endpoints audit actual represented momenta,
support kinetic storage, and per-cell kinetic change/work/heat. Numerical
residuals are reported or refused; they are not deposited as heat. Work and
heat remain separate per-cell outputs. The caller's original caloric owner
must apply total-energy work and retain the associated slip heat once.
Thermal storage, friction calibration, geometry motion, arbitrary 3D support
topology and a new flow clock are outside this mechanical interface.

Workspace allocation occurs at preparation; advancing uses retained buffers
and bounded cancellation checkpoints. A refusal returns no accepted output,
and retry reinitializes all candidates from the supplied original inputs.

Seven focused G0/G1/G3/G4 native tests check unequal-mass/rate analytic
shared motion and individual heat/work, tolerance refinement, symmetric
nonlinear decay, initially comoving cells, Galilean changes and cell order,
clock subdivision, exact single-active parity and workspace reuse, stiff-rate
refusal, and cancellation at every exposed checkpoint with exact retry.
The run compiled original mechanics/drag/controls and fs-math/fs-qty; the
surrounding parent error type was a test fixture. No full crate or workspace
build is claimed.

### Shared branch inertia on the original network source transaction

`NetworkWallOperator::with_shared_branch_supports` binds one optional
`SlidingSupport` per original branch. Every cell on a selected branch reacts
against that one mass and velocity through `SharedRecoilWorkspace`.
Independent per-cell supports may coexist on other branches; overlap is
rejected in either binding order. `None` retains the existing path. Distinct
branch axes are never summed or treated as one spatial body.

Each original source half resolves all selected coefficients at its actual
entering common support velocity. Shared mechanics then supplies per-cell
momentum and mechanical work before the original thermal owner runs with
drag disabled. By default slip loss thermalizes once through total energy minus
kinetic energy. Local wall capacity, original material calorics/conduction, and shared
thermal storage retain their existing ownership. Original transport admission
checks actual relative-slip endpoints, including after shared thermal exchange.
No thermal inventory, mass, mesh volume or flow clock is added.

`with_shared_branch_heat_partition` binds a gas fraction per original branch
cell after common support and original thermal binding. The remainder of each
cell's slip loss is debited from gas before its caloric admission. A finite
local region receives its own share through the original material integral or
native energy coordinate before the existing local heat/conduction owner runs.
An adiabatic local contact instead delivers its share to the bound shared wall:
that owner deposits the compensated sum once, then runs its simultaneous gas
heat solve and existing wall links. Local and shared regions may coexist on a
branch; their stored thermal amounts are already part of its common inertia.
Nonlinear capacity and latent plateaus retain their native source and inventory.
Each heated shared wall belongs to one mechanical support;
other gas ports can exchange heat with that body without claiming its inertia.
Internal energy or explicit negligible expansion is required for material
friction heating. All-gas partition retains the original default path.

The support receipt records actual local deposits and assigned shared transfers.
The shared owner checks its actual deposition and computes wall storage from the complete
pre-friction and post-thermal endpoints. The network combines that physical
change with each gas debit once; represented storage differences remain visible
in residuals instead of being repaired by modifying energy. Source work, both
source halves, wall-link ordering and rollback use the existing transaction.

Accepted support states and `SupportExchange` receipts publish in the same
transaction as gas and thermal walls. Both source halves accumulate actual
friction-only gas impulse/energy; the common support's impulse and kinetic
change use the complete accepted endpoints. Cancellation, exhausted mechanical
work, original source-domain failures and balance refusals preserve accepted
gas, support, thermal states and reports. Preparation allocates the work arrays;
stepping does not allocate. The fixed-axis/frozen-geometry and local estimated
mechanical accuracy limits of the shared solver remain explicit.
Authored heat-fraction rows move into prepared support storage only at binding;
the default all-gas path keeps no extra per-cell fraction vector.

Seven existing shared-support network regressions use the original mixture,
transport, mechanics and thermal implementations. Six further G1/G4 test
groups exercise nonlinear Cv and wall calorimetry against independent integral
equations, unequal heat fractions and films, latent storage through both source
halves of continuing molecular/wave flow, prior accepted-state refusal/replay,
and physical ownership admission. Two regional material inventories connected
by the existing wall-link owner share one support without duplicating its mass.
Local nonlinear and latent regions, mixed local/shared routing, the constant-C
limit, and accepted-state rollback retain their original material sources.
Full Cargo and optional material-feature verification remain separate.

## Two reacting molecular chambers driving one shared piston

`ReactingEuler1d::step_with_shared_piston_and_dissipation` evaluates both
original D/C/A/C/D field compositions inside the same midpoint-speed root.
Pressure and longitudinal viscous traction on both piston faces contribute
their actual signed impulses and work. Each chamber keeps its own molecular
law, species order, chemical references and reactor controls. First molecular
and chemical halves use the old mesh; second halves use actual moved volumes.
The original molecular owner supplies Fourier heat and species-enthalpy
transport, including Soret/Dufour terms when its law provides them.

A reusable private moving-molecular trial also serves the existing single
piston and thermal piston paths, preserving their stage order and arithmetic.
It returns separate per-end pressure/stress impulse and mechanical work,
physical chemical gas-energy change, and finite thermal storage/reservoir
accounting. Conductive end heat never enters mechanical wall-work arrays.

The shared-piston entry requires matching closed constant-area geometry and
clocks. Both fixed outer support impulses remain explicit. Whole-system
admission independently checks actual gas/piston momentum, physical chemical
gas energy, piston kinetic storage and applied-force displacement work. If
either chamber shortens its acoustic interval, both complete compositions
and the mechanical root restart from accepted endpoints. Only the final poll
publishes both fields, meshes, clocks, original molecular geometry bindings
and piston. Existing workspaces are reused through the nonlinear proposals;
no field clone or new constitutive solver is introduced per root proposal.

Four focused regression groups were added for both oriented viscous tractions
against an independent one-cell RK2 oracle, cumulative physical work/energy
and momentum with distinct sources, complete shortened-duration stage replay,
late cancellation with exact retry, and source refusal with both bindings
unchanged. They were parsed by the pinned formatter but were not compiled or
executed in the resource-constrained editing environment.

The shared interface is impermeable and adiabatic. A thermally conducting
shared piston needs one jointly owned finite thermal inventory; two unrelated
end reservoirs must not be substituted for that inventory. Source validity,
explicit parabolic budgets and the existing caloric-interface limitations
remain in force. Symmetric splitting alone establishes no accuracy order,
entropy theorem or benchmarked performance guarantee.

## Two reacting chambers sharing one finite piston

`ReactingEuler1d::step_with_shared_piston` composes both original chamber
mechanisms with the existing two-field, moving-volume pressure/piston solve.
Each source bank retains its species order, caloric references, original
reactor and chemistry controls. The second source halves use the actual
moved cell volumes. A shortened pressure interval restarts both complete
source/pressure/source paths from their accepted endpoints.

Whole-system admission uses both actual chemical-energy changes, the native
piston kinetic change and applied-force displacement work. The two fixed
outer support impulses remain explicit in momentum closure. Independent
source energy gauges are never mixed, and no extra heat-release term is
introduced. Existing trial fields/reactors are reused; both fields, meshes,
clocks and the common piston publish only after the final cancellation poll.

The native regression cases cover independent exothermic mechanisms and
source references, cumulative physical balances, shared-duration replay,
postcompression source refusal and cancellation/retry after prior motion.
They parsed with the pinned formatter but were not compiled or executed in
the constrained editing environment. This entry is adiabatic chemistry plus
Euler transport; molecular transport, finite-wall heat and leakage retain
their separate owners. Symmetric splitting alone establishes no accuracy
order, entropy-stability theorem or combustion-mechanism validity.

### Exact collective modes for common linear shared-support drag

When every active gas amount has the same positive linear acceleration
coefficient and zero quadratic coefficient, the shared recoil workspace uses
the exact collective/deviation solution. Active gas mean/support slip decays
at `a*(1+sum(m_i)/M)`; differences between gas velocities decay at `a`.
Arbitrary gas masses and unequal initial velocities are retained. The mode
uses O(cells) work and one accepted attempt, including stiff/light-support
cases that would otherwise exhaust the generic resolved-step budget.
Disconnected zero-coefficient cells keep their original state and do not
enter the collective mass.

Per-cell physical slip heat and support work use integrated exponential
moments with explicit cancellation-safe weak-coupling evaluation. In
particular, an initially comoving cell retains its induced second-order
impulse and third-order heat. Work remains separate from heat; no kinetic
roundoff residual is relabelled dissipation. Existing actual endpoint
momentum/energy acceptance, compensated accumulators, cancellation and
private-candidate publication remain authoritative. Other shared systems use
the original bounded refined solver, except for the uniform nonlinear groups
described below. The one-active-cell path retains the original exact pair
calculation.

The original network support owner receives this exact mode through its
existing workspace call. No network clock, source splitting, caloric storage
or material-coefficient ownership changes. Source coefficients remain frozen
over the original supplied interval; an exact frozen mechanical solution does
not establish accuracy of a changing material law or the complete multiphysics
composition.

Twenty-one focused actual-source native tests passed: four new collective
regressions, seven shared-support regressions, and ten original drag/pair
regressions. Independent high-precision references cover unequal masses and
initial velocities, individual heat/work, boosts and cell order, both sides
of the moment-evaluation switch, initially comoving cubic heat, sub-ULP mass
contrast, overflowed but fully decayed exponents, disconnected cells, and
cancellation at every exposed checkpoint with exact retry and workspace reuse.
The run compiled the actual controls, mechanics, shared solver and drag with
fs-math/fs-qty; unrelated surrounding error variants were explicit fixtures.
The added stiff common-rate network regression was source-reviewed and
formatter-parsed but not compiled or executed. Full Cargo/network/material
verification is not claimed by this selected-source run.

### Exact shared recoil for uniform nonlinear gas groups

If all active gas cells enter with exactly the same actual velocity and
exactly the same frozen linear/quadratic acceleration coefficients, a
positive quadratic law preserves their common velocity. The shared solver
now advances this invariant group through one original exact reduced-mass
pair map using the compensated sum of active gas masses and the one actual
support inertia. Pure quadratic and combined linear/quadratic laws both
use one accepted interval, independent of the collective relaxation rate.
The existing general common-linear mode still handles unequal velocities.

Each gas cell receives its own mass share of physical slip heat and
integrated work in the initial support frame. The original outer owner
restores work using each cell's actual represented impulse and checks actual
gas/support endpoint momentum and energy. No conservation residual becomes
heat. Disconnected cells do not enter the group and remain unchanged.
Group detection uses exact equality, with no tolerance merging different
states or laws; heterogeneous systems retain the original bounded solver.

The implementation uses the existing workspace with no per-advance vector
allocation. The aggregate unit-density coordinate is private mechanics,
not a new gas cell, volume, caloric state, transport clock or material model.
Original network source owners receive the mode through their unchanged
shared-workspace call. Coefficients remain frozen over the supplied source
interval; material evolution and complete source/flow splitting accuracy
retain their original ownership.

Four new actual-source native regression groups check independent Riccati
slip decay and individual heat/work for both slip signs and pure/combined
laws; unequal density/volume and mass partition; boosts, ordering and clock
partition; stiff/light-support one-interval success and exact heterogeneity
fallback; disconnected cells; and every cancellation checkpoint with exact
retry and retained workspace storage. All 25 selected-source native tests passed, including
the existing common-linear, shared, drag and pair regressions. Full Cargo
and network/material integration execution remain outside this native run.

## Suspended pistons on the original resolved pressure field

`Euler1d::step_with_suspended_piston` and
`MixtureEuler1d::step_with_suspended_piston` add an explicitly authored
`PistonSuspension` to the existing moving-volume RK2 and midpoint-inertia
owners. A stationary axial spring uses nonnegative stiffness, a declared
zero-force world position and the actual represented displacement midpoint.
A viscous damper uses nonnegative damping and the same solved midpoint
velocity. Both loads enter every original pressure/velocity trial; neither
an explicit force kick nor a uniform chamber approximation replaces the
resolved pressure field. The gas mesh remains the sole position owner.
Coefficients and anchor are held over each reported step. Changing stiffness
or rest position between calls can change spring storage outside that step;
the caller owns that constitutive or anchor work.

The accepted report retains the complete native gas/piston report and
separates initial/final spring storage, cancellation-safe quadratic storage
change, viscous loss `c*w²*dt`, anchor impulse and authored drive work.
The original nested additional-force work includes spring/damper work;
`applied_force_work_j` isolates the held external drive. Damper loss leaves
the mechanical system and is not deposited into gas internal energy.
Independent gas-plus-kinetic-plus-spring energy and mechanical-work gates
include the same exterior pressure work and refuse before any publication.
The time-integrated damping loss is independent of represented displacement,
so geometry rounding cannot silently manufacture energy closure.

Every trial reuses the native field/source/geometry workspace. CFL and source
refusals restart the complete attempt at a shortened duration, and a final
cancellation checkpoint precedes joint gas, mesh, clock and velocity commit.
The held-force entrypoints retain their original arithmetic; zero stiffness
and damping reproduce them exactly. Evolving-mixture sources retain each
species and its original formation-energy reference, including exact absences.
These one-sided entrypoints supply a linear stationary suspension, not Coulomb
friction, contact stops, moving anchors, heat return, reacting-piston composition
or an automatic mechanical accuracy guarantee. The requested duration must
resolve the physical spring/damping and gas timescales.

Focused regressions exercise independent midpoint-force solutions, repeated
pressure-driven trajectories with full gas/kinetic/spring/loss and momentum
balance, preloaded equilibrium, exact free-piston reduction, variable-Cp
species and energy-gauge invariance, and early/final cancellation rollback.

## Shared suspended pistons with chemistry and molecular transport

`MixtureEuler1d::step_with_shared_suspended_piston` applies the same original
`PistonSuspension` to both resolved pressure fields around one inertial wall.
Every joint speed proposal retains both moving-volume RK2 fields and their
original independent species axes, caloric sources and energy references.
Both actual face-pressure impulses supply piston motion and conjugate gas work;
there is no imposed exterior-pressure work for an already represented chamber.
Spring force uses the actual shared displacement midpoint, and viscous damping
uses the same solved velocity. Existing retained workspaces own every attempt.

`ReactingEuler1d::step_with_shared_suspended_piston` carries this mechanical
solve through both original C/P/C source compositions. Its
`ReactiveSharedPistonStep<P>` retains the complete chemistry reports; the
suspended mechanics type describes the pressure stage, while the outer
energy/momentum residuals cover the complete chemical step. The default
mechanics type and held-force entrypoint keep their existing report semantics.
`step_with_shared_suspended_piston_and_dissipation` includes both original
D/C/A/C/D chambers inside every speed proposal, with actual pressure and
longitudinal molecular-stress impulse/work from both sides. First source halves
use accepted volumes and second halves use the actual moved volumes.

Complete admission requires both physical chemical-energy changes plus piston
kinetic change, quadratic spring-storage change and independent `c*w²*dt` loss
to equal authored drive displacement work. Each source's formation reference
enters exactly once. Molecular viscous heating stays in the original gas energy;
mechanical damper loss leaves this system. Momentum includes both fixed outer
supports and the suspension anchor's opposite impulse. A shortened duration
restarts both complete chamber compositions. No gas, mesh, clock, piston or
molecular geometry binding publishes before all balances and final cancellation.

These are sealed constant-area axial chambers with one negligible-thickness
piston, a stationary suspension anchor and fixed outer ends. They add no leakage,
piston heat transfer, contact stop, moving anchor or general FSI accuracy claim. The held
coefficient/anchor and caller-owned between-step work conventions remain those
of `PistonSuspension`; the caller resolves the mechanical and gas timescales.
Six focused regressions cover independent preload and midpoint-motion oracles,
full storage/loss/momentum trajectories, distinct source gauges and species,
chemical energy, an independent viscous-gas-heating oracle, exact zero-law
reduction and final cancellation with both molecular bindings reusable.

## Pressure-only outgoing gas boundaries

The opt-in `FrozenGas::pressure_outlet_flux`, `pressure_outlet_normal_flux`
and `pressure_outlet_spatial_flux` APIs impose a nonnegative absolute downstream
pressure on the same original caloric wave owner used by compact reservoirs.
A compressive shock follows the source Hugoniot; an expansion follows the source
entropy and acoustic invariant. The boundary samples the wave at the stationary
face: a departing wave retains the interior state, and a straddling expansion
uses its supported sonic state even when the unused far-field tail is colder
than the source domain.

Exactly zero downstream pressure denotes exterior vacuum for frozen,
evolving-species and native real-gas pressure outlets. The boundary still uses
a positive-pressure, positive-density sonic face on the original caloric/EOS
source, or the unchanged interior trace when all waves already depart.
`FrozenGas::nozzle_exit` likewise accepts zero backpressure and retains the
source-owned stagnation-enthalpy sonic throat. The vacuum branch precedes any
entropy pressure ratio, so it never evaluates a logarithm or an EOS at zero
pressure. If the authored source interval ends before the required sonic
state, evaluation refuses; a cold endpoint cannot replace the sonic state.
No entering vacuum gas state, internal vacuum, exterior expansion solution or
rarefied-flow model is supplied.

The outgoing contact supplies density, temperature and all tangential velocities.
Pressure traction and advected total enthalpy, including every kinetic component,
enter one conservative face flux. The original sensible-energy gauge and
full-trace precision admission remain in force. Primitive fast paths retain
source admission without decoding a freshly reconstructed energy.

A negative contact velocity is refused: pressure alone cannot determine incoming
temperature, entropy or composition. Callers needing reversal must use a complete
reservoir state. These local normal-wave boundaries do not claim nonreflection,
a multidimensional exterior solution, molecular heat/stress boundary data, or
an independently integrated reservoir. Boundary fluxes are oriented along the
supplied coordinate, with `outward_normal` equal to -1 or +1.

Focused G1 regressions compare subsonic compression and expansion, sonic limiting
and supersonic departure against independent gamma-law relations and variable-Cp
entropy/invariant integrals. G3 covers normal reversal, formation-reference
invariance, both shear components and primitive/conserved parity. G4 covers
backflow, invalid data, source precision and cancellation. Vacuum regressions
extend the independent sonic relations, compare nozzle mass flow against the
stagnation choking relation, retain low-positive-backpressure saturation,
exercise unsupported-source refusal and check conservative RK discharge for
frozen and evolving-species fields.

## Heat exchange during driven phase-mixture volume changes

`material::phase_thermal::exchange_driven_mixture_parcel` connects the original
`DrivenMixtureEquilibriumParcel` to the same reservoir, finite constant-capacity
wall or original caloric wall as the isochoric/isobaric thermal entries. The
mechanical caller supplies the final volume and integrated boundary work BY the
mixture. Every heat trial resolves original U,V phase equilibrium with
`Delta U = Q - W_by`, while the shared scalar owner closes
`Q = dt G (Twall(-Q) - Tgas(Q, V_final, W_by))`. No equilibrium heat capacity,
extra latent source or approximated pressure/work trajectory is introduced.

Before choosing the passive heat direction and bounds, the owner solves the
mechanically driven zero-heat endpoint. The initial gas temperature is not a
valid thermal reference for this step: compression can induce cooling when the
wall and gas start at the same temperature, or reverse their original heat-flow
direction. That reference solve consumes the same whole-call heat/PT budgets
and contributes to the returned pressure/EOS counts. It is reused for zero-heat
queries when the numerical controls agree; query counts remain explicit and
PT work is not counted twice. Source brackets must admit the zero-heat driven
endpoint as well as every remaining heat-bracket endpoint.

The accepted report retains supplied heat Q, represented Delta U and prescribed
outward work independently. `mechanical_work_to_gas_j` is `-W_by`, while
`heat_to_gas_j` is the solved Q. Both the source's `first_law_residual_j` and the
outer gas/wall energy residual retain caloric reconstruction defects; neither
turns a numerical defect into a heat source. Zero conductance or zero thermal
duration therefore reports exactly zero heat, but still applies the supplied
volume/work change. Only a complete mechanical and thermal no-op skips phase
work. Mass and original species inventory remain fixed. Gas state, physical
volume and finite-wall storage publish together after balances and final
cancellation; refused source queries, exhausted work and cancellation publish
none of them.

The common wall and paired-mixture first-law checks retain all four independently
represented energy/work terms in a compensated TwoSum reduction. Small heat
cannot disappear when combined with large opposing storage and boundary work.
The physical joule tolerance remains tied to heat/flux, without inflation by
the large mechanical terms. Nonfinite or overflowing reductions refuse;
focused exact-binary accounting regressions check both residual rejection and
correct closure under power-of-two rescaling.

This remains a first-order, spatially uniform, fixed-conductance endpoint
coupling. The caller owns the mechanical trajectory, the actual integral p dV,
its entropy/dissipation accounting, geometry and time refinement. Arbitrary
prescribed volume/work does not establish a realizable piston trajectory or a
combined entropy-production theorem. The original finite phase-screen evidence
and material domains retain their existing scope.

Five focused G0/G1/G3/G4 groups in `phase_equilibrium_thermal` check independent
compression calorimetry with equal initial temperatures, reversal into finite
and latent walls, mechanical evolution at zero thermal exposure, independent
species-datum invariance, whole-call budgets, source refusal and final
cancellation with deterministic retry. They use the existing synthetic binary
constitutive fixture and do not claim experimental-mixture validation.


## Opt-in pressure equilibrium at transported mixture contacts

`MixtureEuler1d::step_with_contact_equilibrium` adds a bounded alternative for
stationary, constant-area inviscid transport. The default `step` remains fully
conservative and retains its unequal-caloric contact-pressure limitation.
The option uses the same species and momentum fluxes, HLLC/HLL wave closure,
canonical species order and original NASA-9 calorics.

A recipient is eligible only when both actual traces on every incident face
match its original pressure and full velocity within the caller's tolerance.
Velocity differences are scaled by sound speed, never by the bulk speed.
HLL fallback, physical boundary cells and later-stage mismatches disable the
correction. Pressure or velocity variations exceeding that tolerance retain the
conservative update.
It reconstructs energy from the transported species and represented momenta
at the original pressure after each forward-Euler stage and the final SSP
average. All source-temperature and kinetic-energy precision gates still apply.

`ContactWorkspace` is prepared once under an explicit byte budget; subsequent
attempts allocate no numerical storage. `ContactEquilibriumControl` requires
a finite positive relative tolerance at most 1e-6 and a finite nonnegative
absolute energy-defect budget per accepted step. The `ContactCorrection`
receipt reports the signed energy defect, sum of absolute corrections and
number of distinct corrected cells. Stage corrections receive weights
1/2, 1/2 and 1, respectively, and each uses the physical cell volume.
Opposing corrections cannot conceal work from the absolute budget.

The native energy residual continues to include this numerical defect.
`corrected_energy_residual_j` subtracts the separately reported signed defect;
the correction is not physical heat or mechanical work. Budget failures halve
the step without forcing MUSCL to first order. Refusal or cancellation leaves
accepted cells and time untouched. Source-inversion hints are private scratch.

This Estimated option does not establish fully discrete entropy stability,
general multicomponent shock accuracy, or equilibrium preservation outside
admitted contacts. Variable-area, body-force, reacting, molecular, piston and
network compositions retain their original transport methods. Moving-wall
transport has the separate opt-in extension below.

## Pressure equilibrium in moving-wall mixture transport

`MixtureEuler1d::step_with_moving_walls_and_contact_equilibrium` extends the
same source-owned contact projection to constant-area closed intervals with
prescribed affine face motion. It retains the original extensive ALE update,
canonical species fluxes, laboratory momentum, moving-wall pressure and
physical `p_star * wall_velocity` work. It uses the same `ContactWorkspace`
and caller-authored absolute energy-defect budget; no additional per-attempt
numerical storage is allocated.

Every actual incident trace must match the recipient's original pressure and
laboratory velocity, with velocity tolerance scaled by the source sound speed.
Mesh speed never enlarges this tolerance. Physical wall cells, resolved
pressure/velocity gradients and HLL fallback stencils remain uncorrected.
The predicted field is projected before its face solve, and the accepted
candidate is projected before source admission and publication. Species,
momentum, geometry and the source's separate species-energy references are
not projected.

For each cell, let `V1` be its moved physical volume and `dE1`, `dEf` the
represented predictor and final energy-density corrections. The implemented
extensive RK2 map satisfies
`Q1 = Q0 + dt*L0 + V1*dE1` and
`Qnext = (Q0 + Q1 + dt*L1)/2 + V1*dEf` on the energy component.
Its signed numerical defect is therefore
`sum(V1*dE1/2 + V1*dEf)`; the absolute budget sums the absolute values of
these same contributions. Both projection weights use moved volumes, including
compression and expansion. The native energy residual retains the defect, and
`corrected_energy_residual_j` removes only this separately reported numerical
change. Boundary pressure work is neither removed nor relabeled as correction.

Budget failures halve the step without lowering the authored reconstruction
order. Each retry resets the contact mask and accounting. Failed or cancelled
attempts leave the accepted species, momenta, energy, complete mesh and clock
unchanged. Zero wall speeds delegate to the stationary contact method,
preserving its field and report exactly; the default moving method remains
fully conservative. Body forces, chemistry, molecular transport, finite piston
inertia and network composition are not enabled by this entry point.

The local uniform-pressure/velocity criterion and the need to expose the
conservation tradeoff follow Abgrall and Karni, *Computations of Compressible
Multifluids*, JCP 169 (2001), section 3.2, DOI 10.1006/jcph.2000.6685.
The volume-weighted accounting above is for this implemented ALE map. This
Estimated extension does not assert general shock accuracy, entropy stability
or pressure preservation outside the explicitly admitted contact stencils.

## Admitted contact pressure driving a finite-inertia piston

`MixtureEuler1d::step_with_piston_and_contact_equilibrium` and
`step_with_suspended_piston_and_contact_equilibrium` compose that same ALE
projection with the existing midpoint piston solve. Every scalar speed proposal
starts from accepted gas/geometry and resets the supplied `ContactWorkspace`,
including bracket endpoints and halved retries. The final evaluated proposal
alone supplies the accepted field, pressure impulse, pressure work, correction,
piston speed and mesh. The native no-correction APIs retain their arithmetic.

The public gas report is `ContactEquilibriumStepReport<StepReport>` inside the
original `PistonStep` or `SuspendedPistonStep`. Its native `transport` retains
actual gas before/after energy and raw conservation residual; its `correction`
retains signed and absolute numerical changes. Both outer
`total_energy_residual_j` accounts subtract only the signed numerical defect
before acceptance. In particular, the admitted free-piston balance is
`Delta E_gas + Delta K_piston + W_exterior - W_drive - Delta E_numerical`.
The suspended account additionally includes the actual spring-storage increment
and nonnegative midpoint viscous loss. Neither pressure work, spring work nor
damper loss is relabeled as a contact correction or gas heating. The absolute
numerical-defect budget and independent mechanical tolerances remain separate.

Species and momentum remain conservative. Physical wall cells, resolved
pressure/velocity gradients and HLL fallbacks remain ineligible for projection;
finite wall inertia does not relax the original contact criterion. Budget
failure halves the coupled interval without forcing MUSCL to first order.
Refusal/cancellation leaves all accepted gas, geometry, clock and piston state
untouched, and the caller's contact scratch is reusable. No per-trial numerical
allocation or field clone is introduced.

Focused G0/G1/G4 regressions cover contact pressure during compression and
expansion, both reconstruction orders, native-map identity on resolved
gradients, final-trial replay against prescribed corrected motion, exact
zero-suspension equivalence, large species-energy-reference invariance,
gas/kinetic/spring/loss/defect closure, budget retries and early/late cancellation.
This remains the same Estimated constant-area, one-piston, nonreacting model.
The shared two-chamber extension has the separate entry points below.
Thermal/molecular/reacting composition, shocks, fully discrete entropy and
general FSI stability are not established here.

## Independent contact pressure on both sides of an inertial partition

`MixtureEuler1d::step_with_shared_piston_and_contact_equilibrium` and
`step_with_shared_suspended_piston_and_contact_equilibrium` run both original
ALE fields with their independently admitted contact policies inside every
joint midpoint-speed proposal. Each chamber keeps its own species/source,
resolution, contact pressure and `ContactWorkspace`. The policies preserve
admitted contacts inside each field; they never equate pressures across the
impermeable partition. Unequal resolved chamber pressures continue to drive
the original piston and optional spring/damper.

`SharedPistonGasReport<R = StepReport>` now carries either the original reports
or one `ContactEquilibriumStepReport<StepReport>` per side. The default APIs
and the conservative helper used by reacting owners retain their signatures
and physical arithmetic. The opt-in reports retain both raw gas-energy changes
and both numerical corrections. Each side enforces its own absolute defect
budget, even if its signed correction opposes the other side's correction.
The outer coupled energy account excludes only their compensated signed sum:
`Delta E_left + Delta E_right + Delta K_piston - W_drive
 - Delta E_numerical_left - Delta E_numerical_right`.
Suspension storage and loss retain their original additional accounts.
Actual interface pressure impulses supply both piston motion and conjugate
gas work. Neither side receives a fictitious exterior-pressure work term.

Both workspaces reset before each paired speed proposal. A refusal in either
field restarts both from their accepted states at the shortened common
duration; a contact-budget retry retains MUSCL. Species, momenta, complete
meshes, equal clocks and piston velocity publish at one final cancellation
boundary. The source/RK/geometry/contact buffers are reused without per-trial
numerical allocation or field cloning. Physical wall cells, resolved gradients
and fallback stencils retain the original conservative map.

Focused G0/G1/G4 checks use different species axes and resolutions, independently
shifted energy datums and unequal chamber pressures. They cover pressure
preservation during motion in both directions, coupled mass/impulse/work/defect
balance, fine-side CFL shortening and replay of both accepted fields, exact
zero-suspension and non-contact limits, spring/damper accounting, independent
budgets with opposing signed defects, and cancellation/refusal of either side
without partial publication. This remains the bounded constant-area sealed
two-chamber Euler model. Thermal/molecular/reacting contact correction, shocks,
fully discrete entropy and general FSI stability remain outside this extension.

## One finite thermal inventory conducting through a shared piston

`SharedPistonThermal` binds the two original chamber end cells to one finite
wall store. `MixtureEuler1d::step_with_shared_piston_and_heat` composes its
actual shared caloric solve with the existing two-field pressure/inertia
operator. `step_with_shared_piston_heat_and_contact_equilibrium` additionally
retains each chamber's independently admitted material contacts and explicit
numerical energy budget. Each face has its own authored complete conductance
`G = h A` [W/K]; the mechanical cross section is not multiplied into it again.
Species do not cross the partition. Thermal mass is counted once and adds no
second mechanical mass to the separately supplied total piston inertia.

The source stages are H(dt/2), A(dt), H(dt/2). Each H uses the original common
`wall_exchange::shared::Group` root to solve both gas heat equations against
ONE wall endpoint. Its calorics borrow each gas's current partial densities,
original NASA source and physical end-cell volume. The first heat stage uses
old volumes, and the second uses the actual moved volumes after pressure work.
The finite wall's source is either an authored constant capacity or its
original material internal-energy chart; an explicit negligible-expansion
solid binding is also admitted. Arbitrary source enthalpy refuses because
this rigid translating wall represents no deformation/pressure-volume work.
Source validity and the original nested caloric work budget remain enforced.
The source-work budget covers both heat halves of each complete attempt and
resets on bounded duration retries; accepted counts exclude discarded attempts.

The first heat stage is independent of the proposed midpoint velocity and is
computed once per duration attempt. Its heated gas fields feed the original
joint pressure root, so heat changes actual piston forces, velocity and work.
Contact targets come from those heated fields, and physical heat gradients
remain unprojected. If pressure transport shortens its duration, both heat
stages and both fields restart from accepted state at that common duration.
Two private flow workspaces and the original heat/RK/geometry scratch are
reused; no full field or caloric table is cloned per velocity proposal.

`SharedPistonThermalStepReport` exposes complete before/after gas totals,
the native middle pressure/work reports, per-chamber heat and the one wall
storage change. Before publication, the enclosing account independently checks
`Delta E_left + Delta E_right + Delta U_wall + Delta K_piston - W_drive
 - Delta E_numerical_left - Delta E_numerical_right`.
The conservative entry has zero numerical defects. Heat never enters the
pressure-work arrays, and the wall's actual represented storage is recomputed
from its accepted and candidate source states. Both gases, complete meshes,
clocks, piston, wall inventory and heat report publish after one final poll.
Refusal or cancellation preserves all accepted physical state and permits retry.

Focused regressions exercise independent three-body backward-Euler equations,
old/moved-volume heat balances, pressure-driven motion from hot and cold walls,
the exact insulated native limit, stiff implicit heat, original nonlinear wall
storage and energy datums, full replay after acoustic shortening, late
cancellation/retry and source-work refusal. Advecting caloric contacts retain
their own pressure while real heating changes pressure at the piston.
A fixed-spatial-discretization time-refinement check follows heat-driven
piston velocity and gas pressure against a finer continuation.
The wall is spatially uniform with constant authored face conductances, and
both outer gas walls are adiabatic. Backward-Euler heat halves establish neither
second order nor a simultaneous conjugate thermomechanical solve. Molecular
transport, chemistry, suspension heat, leakage, thermal expansion, through-wall
temperature gradients, entropy stability and general FSI stability require
their separate owners and are not enabled by these entry points.

## Suspended shared pistons return physical damper heat to one wall

`MixtureEuler1d::step_with_shared_suspended_piston_and_heat` closes the
spring/damper, finite-wall and two resolved pressure fields in one accepted
step. The original `PistonSuspension` enters every existing shared midpoint-
velocity proposal. The caller additionally supplies a finite fraction `f` in
`[0,1]` of the physical viscous loss to return to the one piston thermal
inventory. A damping coefficient alone never selects a thermal destination.

The stages are `H(dt/2), A_with_suspension(dt), deposit, H(dt/2)`. Only an
accepted pressure/mechanical interval pays its actual nonnegative
`D = c*w_mid^2*dt`. The original `Group::after_state` deposits `f*D` through
the wall's nonlinear calorics before the second H; no constant capacity tangent,
second wall mass, additional inertia or empirical gas-heating term is introduced.
Both gas contacts then read that wall temperature and the actual moved volumes.
Returned loss changes subsequent gas pressures and piston motion. The fraction
`D-f*D` leaves the represented thermal system as explicitly reported loss.

The complete suspended account is
`Delta E_left + Delta E_right + Delta U_wall + Delta K + Delta U_spring
 + (D-f*D) - W_drive - Delta E_contact_left - Delta E_contact_right`.
It uses actual endpoint gas and wall storage; returned heat is already inside
the wall/gas account and is not debited again as external loss. The nested
native piston report retains its spring/damper force-work convention and
subtracts `f*D` once from its enclosing gas/wall energy balance. The original
pressure impulses and conjugate gas work contain no conductive or damper heat.
`SharedPistonHeatReport` separately exposes returned heat, its represented
deposition residual and external damping remainder. Its heat-only residual is
`Q_left + Q_right + Delta U_wall - f*D`, not an apparent energy leak.

`step_with_shared_suspended_piston_heat_and_contact_equilibrium` retains the
two independently bounded original contact corrections. Targets come from
the fields after H1; physical heat gradients are never projected away.
Both gas fields, meshes, clocks, piston, wall inventory and last heat report
publish together. Acoustic shortening, failed heat/source admission and
cancellation restore the accepted state before any retry; discarded root
proposals and longer intervals cannot heat the wall. The original heat work
budget includes the caloric deposit. Zero returned heat skips that source
evaluation, and all original pressure/heat workspaces are reused.

Focused G0/G1/G2/G4 regressions cover independent midpoint mechanics and
three-capacity heat equations, actual pressure/motion response to returned
heat, insulated and zero-suspension native limits, nonlinear wall energy and
reference shifts, contact-defect accounting and a duplicate-heat falsifier,
late cancellation after prior motion, acoustic replay, source-work refusal,
and convergence of pressure, motion and wall temperature under time refinement.
The suspension and its anchor remain fixed during a step. Damper heat from
the accepted A does not retrospectively alter that A's force. This split,
spatially uniform-wall model adds no temperature-dependent damping law,
resolved through-wall transport, chemistry, entropy inequality, second-order
thermomechanical accuracy or general FSI stability claim.


### Heat and species diffusion in fixed varying-area channels

`multicomponent::dissipation::DissipationOperator::new_diffusion_only` admits
stationary quasi-one-dimensional meshes with nonuniform widths, varying face
areas and their original integrated cell volumes. The shared RK2 and RKL2
operators use actual area-weighted face divergence, volume-based parabolic
scales and each end's signed area-integrated boundary transfer. Their existing
Maxwell--Stefan mass frame, pressure diffusion, species partial enthalpies,
Fourier heat and authored reciprocal Soret/Dufour coupling remain authoritative.
Uniform-area calls retain their existing arithmetic ordering.

This explicit mode requires zero shear and bulk viscosity at every actual
source evaluation, including later stages; a state-dependent nonzero stress
refuses transactionally. The original full-stress constructor still requires
constant area. This adds no variable-area viscous-stress closure, lateral wall
transport, abrupt-junction mixing model or moving-area geometry. Mesh producers
retain responsibility for physical cell-volume integration. The prepared
workspace is reused without new per-step buffers or successful-path allocation.

Five authored regression groups exercise both integrators against independent
two-volume heat/species relaxation and unequal-end-area thermal-resistance
balances, plus uniform Soret/Dufour parity, late constitutive refusal,
cancellation and workspace reuse. Source tests and isolated arithmetic are not
experimental qualification or a claim that the full workspace suite passed.

## Reacting shared pistons retain one finite thermal store

`ReactingEuler1d::step_with_shared_piston_and_heat` composes both original
chemical mechanisms with the existing `SharedPistonThermal` and two-field
inertial pressure solve. The stages are `C(dt/2), H(dt/2), A(dt), H(dt/2), C(dt/2)`.
First source stages use accepted physical volumes; final source stages use
the actual moved cell volumes. Each chamber retains its own species bank,
caloric path, formation-energy references and fixed-U,V reactor. Formation
energy is converted only by that reactor; neither heat nor pressure work is
inferred from a separate chemistry or constant-gamma law.

Both piston faces share one original finite inventory. The original caloric
wall solve supplies each side's conductive heat and represented storage change.
The first C and H stages change the actual pressures driving the midpoint-speed
root. Final chemistry sees the thermally updated compressed/expanded state.
A shortened inner thermal/mechanical interval restarts BOTH complete chamber
compositions from accepted state at the common duration. Native gas, moving
geometry, source, amount and thermal workspaces are reused.

`step_with_shared_suspended_piston_and_heat` additionally uses the original
held spring/damper. Its authored fraction `f` of accepted loss `D=c*w_mid^2*dt`
enters the SAME wall before H2; the remainder leaves the represented system.
Both final heat contacts and mechanisms therefore see the returned heat.
The complete physical account is
`Delta E_chemical,left + Delta E_chemical,right + Delta U_wall + Delta K
 + Delta U_spring + (D-f*D) - W_drive`.
There is no duplicated thermal mass, formation heat, pressure work or damper
deposit. Conductive and returned heat never enter pressure-impulse/work arrays.
The nested native piston report covers the middle thermal/mechanical stages;
`ReactiveSharedPistonThermalStep` separately exposes complete before/after
inventories, both chemical reports, heat-corrected per-chamber residuals,
physical chemical-energy changes and whole-system energy/momentum residuals.

The thermal owner can retain an unpublished candidate for this composition.
Its accepted storage and heat report publish only after BOTH final chemistry
halves, per-chamber chemical-energy balances, whole-system balances and the
final cancellation checkpoint. Late source refusal, exhausted retries or
cancellation leave both accepted fields, meshes, clocks, piston and wall
unchanged. Every subsequent attempt resets the wall candidate from accepted
storage, so rejected steps cannot retain conductive or returned damper heat.

Focused G1/G4 regressions cover reacting compression with finite-wall energy
closure, accepted-duration replay through actual moved volumes after acoustic
shortening, independent formation gauges, an analytic reversible-isomerization
amount trajectory during heating/motion, independent midpoint spring/damper
mechanics and heat fractions, final-chemistry refusal after wall heating, and
cancellation after prior motion followed by reuse.

This retains the Estimated closed constant-area, one-piston, uniform-wall model,
the original mechanism validity and work budgets, and backward-Euler heat
stages. It adds no molecular transport, material-contact correction, leakage,
thermal expansion, wall temperature gradients, temperature-dependent damping,
general second-order accuracy, entropy inequality or FSI-stability certificate.

### Closed molecular state-entropy admission

`DissipationOperator::advance_for_with_entropy` adds an opt-in entropy and
per-species inventory gate to the original axial SSP-RK2 molecular owner.
Periodic fields, stationary adiabatic reflecting walls, and transmissive ends
with zero molecular heat/species flux are admitted. The existing Maxwell-Stefan,
partial-enthalpy, Fourier, viscous and authored Soret/Dufour fluxes still produce
every candidate. Actual admitted temperatures come from the native total-energy
inversion, including kinetic-to-internal conversion; accepted trajectories are
unchanged when the additional gate passes. The flow clock remains owned by its
composing integrator.

The account uses original source entropy intervals and ideal-mixture partial
pressures. A fixed original temperature and species-density basis removes
arbitrary linear species reference terms from the admission quantity.
Double-double density differences and a cancellation-safe relative-entropy
remainder retain weak mixing and finite absent-species limits. The report
separates the admitted thermal/mixing increment, the original source-state
increment, and the reference contribution from numerical inventory residuals.
Each species' continuity residual is independently bounded relative to its
original mass; an originally absent species must stay absent. Conserving only
total mass or elements cannot substitute for this species check.

Every actual candidate must satisfy the caller's absolute entropy allowance.
Entropy or inventory failure uses the existing bounded step-halving policy,
without clipping species, temperature or energy. After all internal steps,
the direct original-to-final entropy and species account must satisfy the same
single allowance before any caller field is published. Internal accepted
substeps therefore cannot accumulate a larger allowance. Cancellation, work
exhaustion and final refusal preserve the complete original field and clock.
Only the opt-in lane allocates entropy temperature storage; it reuses native
candidate inversions and adds no spatial constitutive evaluations.

This is a measured finite-step state-entropy admission, not a decomposition into
constitutive physical production or a fully discrete entropy/invariant-domain
theorem. Spatial and temporal numerical effects remain in its signed result.
In particular, the existing arithmetic partial-enthalpy face can have negative
entropy rate even while conserving physical energy; reducing the step need not
repair that spatial defect. Reservoir chemical/thermal entropy exchange,
moving walls, reactive or advective composition, RKL2 and network/spatial
owners are outside this entry point.

Five focused G0/G3/G4 regression definitions cover native trajectory parity
under all three admitted boundary types and large formation-reference shifts;
independent thermal/mixing entropy with viscous conversion; an energy-conserving
negative-entropy face that must refuse; complete-interval allowance and rollback;
weak mixing, absent species, late cancellation and reuse.

### Entropy admission with stationary molecular reservoirs

`advance_with_reservoirs_and_entropy_for` applies the original axial SSP-RK2
molecular update against the existing source-bearing stationary end baths,
then admits actual state entropy minus actual bath entropy input. Each end
retains its own temperature, species composition and original caloric source.
The native species and energy transfers use the original two stage weights,
physical end areas and half-cell boundary distances. Heat, species enthalpy
and source chemical entropy are counted once; no advective flux is introduced.

For signed transfers into the gas, the bath contribution is
E_ref/T_b - sum_i[(h_i,ref/T_b-s_i,b)*dm_i/M_i]. This equals the native reduced
heat divided by T_b plus the bath partial-species entropy transfer. Bath
enthalpies are those already used by the native energy flux. The stationary
bath has no mechanical-work contribution. Fourier and authored reciprocal
Dufour heat enter through that same energy transfer, with no second heat source.

The state and boundary accounts share one fixed original temperature/density
basis. Original molar internal-energy and entropy intervals cancel arbitrary
formation and entropy references before reduction. The report separately
retains full-source state entropy, left/right full-source entropy input,
inventory-corrected state-minus-input generation and the source-reference
contribution from each species' continuity residual. These are finite-step
state accounts, including spatial and temporal numerical effects, not isolated
constitutive production integrals.

Each species' residual is its actual mass change minus native boundary
transfer. Its relative allowance scales against original mass plus accepted
absolute stage/end throughput; opposing ends or stage reversals cannot cancel
that scale. An initially absent species may enter from a declared bath.
A bath species with zero activity must have zero flux in both actual stages;
a finite chemical potential is never fabricated by flooring a concentration.
A bath with absent species can still exchange heat when those species fluxes
are exactly zero.

The same absolute entropy allowance applies to every native substep and to
the direct original-to-final state-minus-input account. Failed candidates keep
the original bounded halving policy. All final entropy/report arithmetic is
completed before the caller's field is published, including the full-source
reference contribution. Cancellation, work limits, zero-activity refusal and
final entropy refusal leave the complete original field and clock unchanged.
Closed entry points retain their existing semantics. Moving baths, advection,
chemistry, RKL2, and network/spatial entropy remain outside this entry point.

Four focused native regression definitions cover unequal-temperature steady
Fourier baths; composition baths that reduce body entropy while total generation
is positive; original trajectory and formation-datum invariance; admission of
initially absent incoming species; zero-activity refusal and heat-only admission;
and cancellation before final publication followed by workspace reuse.

### Entropy-conjugate axial molecular transport

`DissipationOperator::advance_entropy_consistent_for` selects a source-conjugate
interior face and retains the actual SSP-RK2 state-entropy admission. It uses
the same original gas, Maxwell-Stefan solver, material/kinetic transport law,
geometry and native publication transaction. Compile-time selection in the
shared native loop retains the original faces for existing entry points;
no parallel time integrator, altered public controls or extra owner state is
introduced.

Each species transports its original divided molar enthalpy
h_i*=Delta(h_i/T-s_i^0)/Delta(1/T), with the existing reference-energy convention.
The implementation forms this mean from source caloric/entropy intervals
around the sampled temperature, before large formation references can cancel.
Equal temperatures use the original referenced enthalpy directly. Variable
heat capacities, piecewise caloric intervals and explicit source continuity
remain those of the original thermochemical owner; no constant-Cp or
constant-gamma surrogate is introduced.

The conjugate chemical effort is Delta(mu_i/T)-h_i*Delta(1/T). Its actual
finite-precision remainder is retained. Explicit barycentric balancing
supplies the original Maxwell-Stefan driving-force solve. Authored Soret heats
add their reciprocal-temperature force and the same coefficients multiply
actual species flux in Dufour heat. Arithmetic velocity in viscous work closes
the paired axial stress contribution, including unequal cell widths.

The paired semidiscrete entropy rate consists of Maxwell-Stefan pair friction,
Fourier heat production and axial viscous production. Original row/mass
residuals and floating-point reductions remain limitations of that identity.
The complete native finite-step entropy and species-inventory gates still
apply before publication. This supplies neither an unconditional timestep
bound nor a fully discrete entropy/invariant-domain theorem.

Constitutive T, p and mole fractions retain the original width-weighted
sampling, preserving exact equal traces. The coefficient law is evaluated
once per face, with the original separate query for authored cross
coefficients. Reusable Maxwell-Stefan storage and bounded stack arrays are
retained. The engineering diffusivity estimate includes trace composition
and temperature contrasts; it is not a nonlinear spectral or positivity
certificate. Equal-temperature faces skip unnecessary caloric interval queries.

The admitted interior domain has common species support: each species is
positive at both traces or absent at both. One-sided zeros refuse before a
transport-source query. There is no concentration floor, invented trace mass,
post-solve flux clipping or fallback to the engineering face. Periodic and
stationary zero-heat/species ends retain their native closures. Reservoir and
moving-wall evolution, network/spatial faces, reactions and advective transport
remain separate entry points.

Focused regression definitions cover the original energy-conserving
negative-entropy counterexample and its repaired face and full SSP-RK2 step;
variable and piecewise source calorics; equal/adjacent temperatures; formation-
reference shifts; reciprocal thermal diffusion and viscous work; common
absence; and transactional domain/cancellation refusal. The independent
two-cell initial entropy rate changes from -0.0397353103431 to
+0.0675891418055 W/K without reversing mixing or injecting energy.


### Source-conjugate stationary molecular baths

`advance_entropy_consistent_with_reservoirs_for` evolves the original axial
SSP-RK2 field with conjugate faces at BOTH interior interfaces and stationary
source-bearing reservoirs. The same divided original caloric enthalpies,
logarithmic chemical forces, barycentric Maxwell-Stefan solve, Fourier flux
and authored reciprocal Soret/Dufour terms supply one conservative energy
flux. Each bath layer has the original physical half-cell distance and actual
end area. Arithmetic velocity in stress work gives the paired viscous
production, including a nonzero evolving cell velocity next to a stationary
bath. Bath traces are fixed sources, not extra evolving cells.

The entropy input of each bath is evaluated from the ACTUAL conservative
transfer as E_ref/T_b - sum_i(mu_i,ref/T_b * dm_i/M_i), on the existing fixed
source-reference basis. Relative to bath partial enthalpies, its heat is
E_ref - sum_i(h_i,ref(T_b) * dm_i/M_i). That heat contains the divided-enthalpy
layer correction and the bath share of viscous conversion. The interior
Fourier/Dufour reduced heat alone is not the bath heat; substituting it would
omit part of the actual energy transfer. No additional heat or traction work
is deposited. A stationary bath does no external mechanical work.

Native candidate positivity, per-species continuity, substep entropy and the
single original-to-final entropy allowance all remain in force. The complete
field and final entropy report publish only after those checks. Both closed
and legacy reservoir entry points retain their original face choices through
compile-time selection in the same owner and SSP-RK2 loop. Constitutive
sampling, reusable diffusion workspace and bounded stack arrays are shared;
accounting introduces no second face/source evaluation or per-face allocation.

Every species has common support at each admitted interior AND bath face:
positive on both traces or absent on both. One-sided zero activities refuse,
without clipping, invented trace mass or a switch to the engineering flux.
The legacy reservoir method remains available for its own separately admitted
domain. Moving baths, advective/reactive entropy, nonideal mixtures and
network/spatial owners remain outside this entry point. The semidiscrete pair
identity and measured finite-step gate do not imply an unconditional timestep,
invariant-domain theorem or experimental validation.

Three focused native regression definitions cover a bath-enthalpy
negative-entropy counterexample and its repaired full SSP-RK2 step; both end
orientations; variable heat capacity, unequal molar masses, reciprocal thermal
diffusion, Fourier heat and viscous conversion; formation-reference shifts;
one-sided-zero refusal, late cancellation and workspace reuse.
For the declared hot-bath/cold-cell fixture, the initial entropy rate changes
from -0.1758153103431 to +0.0675891418055 W/K. At dt=1e-4 s the body entropy
decreases while the source-aware total generation is positive, so a closed-body
entropy test would be the wrong admission criterion.


### Source-conjugate RKL2 molecular evolution

`advance_stabilized_entropy_consistent_for` combines the native RKL2
super-time-stepping recurrence with the original-source conjugate molecular
faces. The same original caloric models, barycentric Maxwell-Stefan solve,
Fourier heat, longitudinal stress and authored reciprocal Soret/Dufour terms
remain responsible for the trajectory. Compile-time face selection in the
existing native recurrence preserves the original stabilized entry points;
there is no parallel time integrator or added operator state.

Every native intermediate and final stage retains the original species,
caloric, source-domain and parabolic-scale admission. Actual state entropy
and each species inventory are additionally checked at every physical
macrostep endpoint. Internal Legendre stages are algebraic auxiliaries,
not sequential physical timesteps; stage-by-stage entropy monotonicity is
not asserted. A failed endpoint uses the existing bounded halving retry.

The endpoint account reuses temperatures already inverted by the final native
RHS admission. There is no additional EOS inversion, constitutive evaluation,
per-stage allocation or additional persistent temperature buffer. Candidate
temperatures replace the entropy anchor only after admission, so a rejected
endpoint cannot contaminate a retry. The direct original-to-final entropy and
species account must satisfy the same single allowance before publication.
Cancellation, exhausted work or final refusal leaves the complete field and
clock unchanged; no accepted prefix is published.

This remains an Estimated nonlinear transport lane. The RKL2 negative-real
stability interval does not establish nonlinear positivity, entropy stability
or temporal accuracy for the coupled system. The caller retains the physical
macrostep accuracy control and bounded stage/RHS budgets. Admitted species
have common support at each face, and boundaries are periodic or stationary
with zero heat/species exchange. Reservoir, moving, advective and reactive
evolution remain separate entry points.

Four focused regression definitions cover independent Legendre heat spectra
and actual source entropy, actual RHS counts against SSP, the repaired
unequal-caloric face and its two-stage Butcher trajectory, source-reference
shifts and common absence, and late cancellation/work refusal with reuse.
A positive smooth conductivity fixture stays within native caloric and
parabolic limits yet yields a negative-entropy long RKL2 endpoint; the new
gate rejects it, retains the correct retry anchor, and accepts shorter
physical macrosteps. This is a concrete counterexample to inferring nonlinear
entropy admission from the linear stability interval alone.


### Accelerated source-conjugate transport with stationary baths

`advance_stabilized_entropy_consistent_with_reservoirs_for` connects the
source-conjugate reservoir faces to the SAME native RKL2 recurrence. Both
thermal and species exchange now use the accelerated operator with original
calorics, material/kinetic transport, Soret/Dufour coupling and viscous work.
The bath states are fixed over the requested interval; their entropy efforts
are formed once from the original source and the existing reference basis.

The entropy input and per-species continuity checks use the final signed
left/right energy and species transfers retained by the native conservative
RKL2 recurrence. They are not replaced by SSP weights, a separate quadrature
reconstruction, or reduced Fourier/Dufour heat. Stationary bath entropy is
linear in these actual transfers, so no repeated bath-source evaluation is
required. The extra final RHS admission has zero quadrature weight.

RKL2 has signed effective RHS weights: its two-stage form has weights
[-1/2, 3/2]. `Rkl2Scheme::rhs_weights` obtains the native tableau weights by
reverse accumulation in O(s) work using caller storage. No normalization or
positivity repair is applied. Inventory tolerances use original species mass
plus accepted absolute weighted stage/end throughput. Taking the absolute
value only after stages or ends cancel would understate this numerical scale;
taking absolute values at every recurrence level would inflate it. Neither
shortcut is used. These are numerical quadrature throughputs, not separate
physical evolution of the auxiliary Legendre states.

Each retry clears its trial boundary account. Only an accepted macrostep adds
to the cumulative transfer and throughput. Each actual endpoint and the
direct original-to-final state-minus-bath account must satisfy the same
absolute entropy allowance, while each species retains its independent
continuity bound. Complete source-entropy report arithmetic precedes native
field publication. Late cancellation, source/domain refusal, exhausted work,
or final entropy refusal leaves the caller field and clock unchanged.

This adds bounded O(s) weight work and bounded stack storage, with no extra
RHS, EOS inversion, or per-stage allocation. Original RKL2 work budgets,
physical macrostep accuracy controls and intermediate state admission remain
active. The admitted closure is fixed axial geometry and stationary
source-bearing ideal-mixture baths, with common support on both traces of
every active face. It does not establish nonlinear unconditional stability,
a general positivity theorem, moving/time-dependent baths, or advective,
reactive, network or spatial entropy evolution.

Three focused native regression definitions cover a multistage Dirichlet heat
transient with separate actual bath temperatures and an independent Legendre
mode oracle; signed two-stage binary species exchange, source-reference
invariance and direct accepted-throughput inspection; and late cancellation,
incompatible bath support and reuse. Both thermal and composition fixtures
have negative body entropy change but positive state-minus-bath generation.
The binary fixture has near-zero net species transfer across the complete
boundary while retaining nonzero absolute contributions from both stages
and both ends.


### Source-derived inviscid entropy-conjugate mixture faces

`MixtureGas::entropy_conservative_flux` and `entropy_stable_flux` add original-
source axial faces for the existing conserved species, momentum and referenced
total energy. Every species retains its own molar mass and temperature-dependent,
possibly piecewise NASA-9 calorics. The original EOS supplies both endpoint
acoustic speeds. Reattaching the original energy datum to each species flux
recovers physical total-energy transport, including the dissipative correction.

For mathematical entropy eta=-S, beta=1/T and r_i=R/M_i, the entropy variables
are v_i=(h_i,ref/T-s_i,partial)/M_i-beta*u^2/2, v_m=beta*u and v_E=-beta.
The entropy potential is psi=u*sum_i(r_i*rho_i). The original molar function
phi_i=beta*u_i,ref-s_i^0+R*ln(T) satisfies d(phi_i)/d(beta)=u_i,ref.
Its divided energy e_i#=Delta(phi_i)/(M_i*Delta beta), logarithmic species
densities and arithmetic velocity produce the conservative face:
F_i=rho_i^log*ubar; p#=sum_i(r_i*rhobar_i)/betabar;
F_m=ubar*sum_i(F_i)+p#; and
F_E=sum_i(e_i#*F_i)+u_left*u_right*sum_i(F_i)/2+ubar*p#.
Product jump identities give Delta(v).F_ec=Delta(psi) in exact arithmetic.

The divided energy is evaluated through original internal-energy and entropy
intervals around one sampled temperature. Formation-energy terms cancel inside
those source intervals before rounding; caloric seams retain the source owner's
continuity contract. Equal temperatures use the original referenced internal
energy directly and skip three interval queries. Stable logarithmic ratios,
temperature differences before inverse-temperature subtraction and double-double
reductions retain weak jumps. Exactly identical admitted traces return the
native physical flux of the original rounded conservative state.

The dissipative choice subtracts alpha*Delta(q)/2 in every conserved component,
with alpha=max(|u_left|+c_left, |u_right|+c_right). Convexity of -S on the admitted
common-positive caloric/species domain gives the semidiscrete numerical entropy
contribution alpha*Delta(v).Delta(q)/2 >= 0. This is numerical dissipation, not
a second molecular heat source. `EntropyFaceReport` separately returns this
contribution, the signed conservative identity defect, and the contraction
against the actual final rounded flux. No residual or negative sign is clamped,
and no flux repair is applied to force an entropy receipt.

The successful face path allocates nothing and polls cancellation while walking
the bounded original species axis. The native private interface reuses decoded
endpoint states and compile-time removal of diagnostic contractions; it adds no
EOS inversion. Every species must be positive on both traces or absent on both.
Common absence carries exactly zero flux; one-sided absence has a named refusal,
without mass floors or a switch to another face.

These are source-conjugate semidiscrete faces, not a finite-step invariant-domain
theorem, reconstruction certificate, or guarantee of pressure equilibrium at
moving material contacts. Accepted native time evolution needs its own actual
state/inventory/entropy admission. Viscous stress, molecular heat/species
transport, reactions, moving boundaries and spatial/network geometry remain
with their distinct original owners. The calorically perfect mixture construction
and its positivity/contact limitations are described by Gouasmi, Duraisamy and
Murman, https://arxiv.org/abs/1904.00972; the variable-caloric source-interval
divided energy above is derived for this owner's original EOS.

Four focused regression definitions use independent polynomial entropy
derivatives and source integrals, unequal molar masses, a piecewise caloric seam,
formation gauges of 1e18, physical-energy reattachment, exact consistency,
zero velocity, adjacent temperatures/densities, common absence and cancellation.
For the declared two-species variable-Cp fixture the conservative flux is
[0.22939772576023482, 0.1796030100332838, 2.5951614422181946, 4.539804789164915]
and numerical entropy production is 11.711975129142397 W/(m^2 K).
Source regression definitions and isolated arithmetic are not a claim that the
native Rust tests ran or that this model has experimental qualification.


### Native inviscid entropy-admitted SSP-RK2 evolution

`MixtureEuler1d::step_entropy_consistent` connects the original-source
dissipative inviscid face to the existing axial SSP-RK2 owner. It advances one
accepted native step, returning the original `StepReport` and the shared
`MixtureEntropyReport`; `MixtureEntropyControl` aliases the original state
entropy controls. The requested duration is a maximum: native CFL and bounded
admissibility/entropy retries can shorten the actual returned duration.

This entry point admits first-order space, fixed constant-area geometry, and
periodic or two stationary reflecting boundaries. The native reflected-wave
wall closure supplies wall momentum and no species, energy or external entropy.
Every active species must have common positive support across each interior
face, or be absent on both traces. Open advective ends, reconstruction, moving
or varying-area geometry, forcing, reactions and coupled molecular composition
require their own complete treatment and are not enabled by this entry point.

The original integrator still owns both conservative stages, intermediate
caloric/positivity admission, actual face-speed CFL checks, bounded retries,
extensive balances and final state/clock publication. A compile-time interior
face choice reuses that integrator and the already decoded original states.
Existing engineering, forcing, moving and contact-correction paths retain their
original dispatch. No second integrator or replacement caloric source is added.

The final native candidate loop passes its already admitted temperature and
actual conserved state into the same source-interval entropy account used by
molecular SSP-RK2 and RKL2. The account checks every species' inventory and the
original-to-final thermal/mixing entropy before native success. Entropy or
inventory refusal enters the original halving loop without state clipping,
mass seeding or an engineering-flux fallback. Each retry resets its anchor
temperatures from the unchanged original field. No second candidate EOS solve
or repeated final entropy-interval evaluation is introduced.

A common pure report reducer completes all finite arithmetic before accepting
the account. The original final balance checks and cancellation poll still
precede the field swap and clock advance. Failed or cancelled work cannot
publish a partial field. The entropy allowance applies to this complete
accepted step; it is not a cumulative allowance across separate calls.
Source-reference effects of numerical species residuals are reported separately
and cannot replace the individual continuity checks. Molecular public errors
keep their existing variants through explicit conversion from the shared native
state-account errors.

Four native regression definitions cover an independently computed two-cell
SSP-RK2 pressure/density mode, variable-Cp unequal-mass reflected pressure waves,
large formation gauges and common absence, exact uniform motion, and late
cancellation/open-boundary refusal with successful reuse. The two-cell fixture
has an independent entropy increment of 0.027299800511537853 J/K for its
0.02-second accepted step. These are regression definitions and independent
arithmetic evidence; native Rust execution and typechecking were unavailable.
The semidiscrete face identity and these runtime gates do not establish a
general finite-step invariant-domain theorem or moving-contact pressure
equilibrium guarantee.


### Closed entropy-consistent reacting-flow composition

`ReactingEuler1d::step_entropy_consistent` and
`step_entropy_consistent_with_dissipation` connect the original source-derived
inviscid and molecular fluxes to the existing coupled C-A-C / D-C-A-C-D owner.
The original fixed-U,V chemistry retains its source identity, kinetic laws,
selected implicit method, positivity, element/charge/mass, energy, local-error,
entropy and work admission. Chemical heat release enters through the original
species-energy gauge update exactly once. No second reactor or split integrator
is introduced.

Every accepted component uses one actual duration. If the middle inviscid step
shortens for acoustic or entropy admission, the existing retry loop restores
accepted state and replays all chemical and molecular halves. Both optional
SSP-RK2 molecular halves use original source-conjugate Maxwell-Stefan, Fourier,
viscous and authored reciprocal Soret/Dufour transport. Their finite-duration
entropy gates remain active. Only the middle transport advances physical time.

The complete candidate has a separate original-to-final entropy and species
source check inside that same retryable attempt, before field/clock publication.
Write the existing state basis as S=F+sum(alpha_i*m_i/M_i). Chemical species
changes Q_i are physical, so the admitted whole change is
Delta F+sum(alpha_i*Q_i/M_i). The original source change is
Delta F+sum(alpha_i*Delta m_i/M_i); their difference
sum(alpha_i*(Delta m_i-Q_i)/M_i) is only the numerical continuity-reference
defect. Removing the entire species reference term would incorrectly remove
genuine reaction entropy. Adding chemistry-only entropy to the complete endpoint
would instead count chemical change twice; neither shortcut is used.

Q_i and its absolute throughput are observed from each actual original
before/after chemical cell row, with double-double differences before volume
multiplication and accumulation. Each species residual Delta m_i-Q_i is bounded
relative to original inventory plus the sum of absolute chemical changes across
cells and halves. Production and consumption cannot cancel this scale. A species
absent initially can be created; zero original inventory and zero throughput
still require exactly zero drift. Each retry clears trial transfers and restores
immutable original temperatures. The final state decode also checks the original
chemistry temperature interval after molecular transport, sharing that pass with
the whole entropy calculation.

`ReactiveEntropyStepReport` retains the original conservative coupled report,
direct `ReactiveEntropyReport`, middle flow entropy and optional
`SplitMolecularEntropyReport`. Component receipts do not replace complete-state
admission. One caller allowance applies directly to the entire accepted step;
component allowances are not added into a larger final allowance. Every report
reduction, original balance check and final cancellation poll precedes native
publication. Existing non-entropy compositions retain their original entry
signatures and dispatch.

The new scope is first-order space, fixed constant area, periodic or two
stationary reflecting ends, with common species support at active transport
faces. Open boundary entropy, moving geometry and higher-order reconstruction
are not inferred. Symmetric splitting does not change the selected chemistry
method's order or establish an unconditional positivity/entropy theorem.

Four regression definitions cover reaction-created species and heat-release
duration replay, variable-Cp reacting/molecular flow with a nonzero original
species entropy datum, independent endpoint entropy and energy balances, late
cancellation/reuse and unsupported spatial/boundary scope. An independent
fixed-energy polynomial example gives Delta S=0.09976352648220975 J/K, including
0.04494595720774408 J/K of real chemical source-reference contribution. These
arithmetic checks and source parsing do not claim native Rust execution or
typechecking, which were unavailable.


### Source-conjugate RKL2 in complete reacting entropy evolution

`ReactingEuler1d::step_entropy_consistent_with_stabilized_dissipation` selects
the original RKL2 recurrence for both molecular halves of the same closed
D-C-A-C-D composition. It reuses the preceding explicit entry's scope admission,
chemistry, source-aware whole-state account, native duration restarts and atomic
publication. The inviscid middle step and both molecular halves use their
original source-derived entropy fluxes; no duplicate constitutive operator,
reactor or composition loop is added.

The original molecular RKL2 owner now has a counter-preserving conjugate entry.
One actual RHS counter is retained across BOTH halves and ALL discarded coupled
attempts. Acoustic shortening, entropy rejection or source refusal cannot reset
the allowance. Original stage/end-state calls consume that counter before
evaluation; work exhaustion remains terminal in both molecular and coupled
retry loops. Each retained half report describes its accepted molecular call,
while `StabilizedReactiveEntropyStepReport::operator_evaluations` includes
abandoned coupled work as well.

Every actual molecular macrostep endpoint and complete half interval retains
the original caloric, positivity, constitutive, species and entropy admission.
Internal Legendre stages remain algebraic auxiliary states; they are not
misclassified as consecutive physical entropy-producing time steps. The whole
reacting state is still independently admitted after the complete composition.
Caller-selected macrostep and chemistry controls remain necessary; no global
temporal-error, nonlinear positivity, contact-equilibrium or combustion theorem
is inferred from symmetric splitting or the RKL2 linear stability interval.

Two focused regression definitions compare a diffusion-dominated smooth coupled
mode with the explicit source-conjugate split, check independent endpoint
entropy and physical-energy balances, exercise extended stage counts and actual
RHS work, and exhaust a three-call budget after an abandoned heat-release
attempt. That refusal must leave the public field/clock unchanged and the same
owners reusable. The tests are definitions, not executed performance results:
source/syntax review completed, but native Rust typechecking and test execution
were unavailable.


### Complete-duration entropy-consistent reacting evolution

`ReactingEuler1d::advance_stabilized_entropy_consistent_for` and
`advance_stabilized_entropy_consistent_to` connect the preceding source-conjugate
reacting step to the existing full-duration RKL2 interval owner. The original
loop, compensated physical ledger, exact target handling, reused conserved
rollback buffer and aggregate work controls remain authoritative. No duplicate
interval integrator or complete flow clone is introduced.

One original source-entropy account anchors the entire requested window.
Each native coupled step retains its own unchanged retry anchor and all
component/complete-step checks. Only after native success and interval clock
validation does its exact double-double chemical species change and absolute
per-cell/per-half throughput merge into the window account. Discarded attempts
contribute work, but no chemical source, to the accepted window. These transfers
are independent of each step's entropy basis; rounded per-step entropy reports
are not summed to reconstruct the final state.

At completion the same original source-interval reducer compares the actual
final field with `interval_initial`, using the immutable original temperatures
and source reference basis. It restores genuine chemical source entropy,
separates only numerical species-source reference defects and admits each
species against original mass plus accepted absolute chemical throughput.
ONE original-to-final entropy allowance applies to the complete window,
in addition to every step's local admission. Local allowances do not accumulate
into a larger final allowance. Endpoint scratch resets preserve merged accepted
chemical transfers and the original anchor.

Both molecular halves, every RKL2 stage, all discarded coupled attempts and all
accepted physical steps share the existing interval RHS/coupling budgets.
Each terminal remainder is actually integrated. A duration unrepresentable at
the current clock refuses; the target-time companion lets the caller explicitly
select a represented endpoint. No clock snap or epsilon remainder is substituted.
A no-op still admits the original operator, source and closed entropy scope.

The complete entropy gate runs inside the original rollback-protected closure,
after the physical interval report is constructed and before the final
cancellation poll. Source or entropy refusal, exhausted work, stagnant time or
late cancellation restores the ENTIRE original accepted field and clock,
including any internally completed prefix. The same owners remain reusable.

`StabilizedReactiveEntropyIntervalReport` returns the original physical interval
report, direct complete-window `ReactiveEntropyReport`, smallest accepted
coupled-step entropy change and largest accepted-step species residual. The
complete endpoint also passes its own direct species check. The original
first-order, fixed constant-area, closed-boundary, common-support limitations
remain; no general temporal-error, positivity or combustion theorem is inferred.

Focused regressions cover multiple-step variable-Cp reacting flow with a
nonzero species entropy datum, independent complete-state entropy/source
balances, global work exhaustion after an accepted prefix, late cancellation
and reuse, no-op/represented-time behavior, and the unchanged whole allowance.
The latter uses actual source states with individually allowable entropy
changes of about -6e-6 J/K each and a direct total of -1.2e-5 J/K: a single
1e-5 J/K whole allowance must refuse that accumulated decrease. Syntax/source
review and independent arithmetic were performed. Native Cargo tests,
typechecking and the DSR gate were not run; disk pressure and unavailable
remote build runners prevented a supported native verification run.

## Photon-conserving atomic gas evolution

`chemistry::photoionization::advance_monochromatic_photoionization` evolves
a rigid, uniformly mixed column under a constant monochromatic photon number
flux. Its original atomic source supplies the event cross section and chemical
gap. The actual neutral inventory supplies Beer-Lambert opacity. The update is
analytic: `expm1(tau1) = expm1(tau0) * exp(-sigma*flux*duration)`.
There is no reaction time-step or nonlinear-iteration loop in this owner.

Stable thin/thick branches calculate the smaller absorbed or surviving
population directly, preserving tiny ionization increments and representable
neutral tails without overflowing at large optical depth. Each new ion and
electron corresponds to one absorbed photon. Transmitted photons supply no
heat. Original chemical storage and the changing electron sensible pool close
against absorbed photon energy within an explicit numerical joule allowance.
Conserved formation-reference offsets cancel before the physical energy
difference. Heavy temperature and all spectator amounts remain unchanged.

Exact zero initial electrons are allowed: photoionization supplies their
physical source. The endpoint is the existing `FrozenIonizedGasStateV1`, with
its original thermodynamics, amounts, separate temperatures and charge checks.
Unsupported electron temperatures, source domains, arithmetic and cancellation
refuse before publication; no temperature floor or artificial seed is added.

This is a perfectly mixed, quasi-static monochromatic column model. It does
not resolve a spatial ionization front, light-travel delay, recoil, secondary
ionization, excited states, radiative recombination or a nonthermal electron
spectrum. Immediate photoelectron thermalization is an explicit assumption.
Radiation entropy is not modeled, so closed-material entropy inequalities are
not applied to this driven subsystem. Discrete spectral source evaluation
exists in L1; different photon cross sections cannot be collapsed into this
single-group exact transient.

Focused regressions specify zero-seed creation, original-source energy and
mass/charge accounts, photon-limited thick columns, thin exponential depletion,
tiny surviving populations, exact interval composition, conserved-reference
invariance, source refusal and cancellation. The scalar kernel includes
independent high-precision optical-depth reference values.

## Layered photoionization and evolving opacity

`chemistry::photoionization::column::advance_monochromatic_column` transports
one fixed-area monochromatic beam through ordered, individually mixed gas
cells. Each cell retains its own initial state and original atomic source.
The upstream analytic solution supplies escaped integrated photon fluence to
the next cell. Since the admitted photoionization endpoint depends only on
cumulative incident fluence, one analytic solve per cell accounts for evolving
upstream opacity without freezing it or introducing chemical time substeps.
Reported downstream fluxes are fluence-equivalent means, not instantaneous
illumination histories.

The single-cell owner evaluates outgoing photons directly with stable
thin/thick expressions. It retains representable transmitted populations that
would disappear in incident-minus-absorbed subtraction, including scaled tails
for which unit transmittance underflows. Those photons can create downstream
electrons and sensible energy through the original atomic source.

An explicit maximum cell count bounds admission and allocation. Compensated
incident-to-exit reductions check the whole-column photon account and a single
unchanged absolute energy allowance. Indexed source failures, arithmetic
refusals and late cancellation publish no partial endpoints; caller states
remain immutable. Positive unrepresentable downstream mean flux refuses.

The exact serial endpoint statement requires fixed geometry, fixed source
cross sections, quasi-static propagation and the parent owner's immediate
photoelectron thermalization. Recombination, cooling, collisional chemistry,
scattering, other absorbers, motion and light-travel delays are absent during
the interval; structure within each cell remains unresolved. This is not a
general radiation-hydrodynamics or multi-frequency ionization-front solver.
Controls include independent homogeneous-prefix and heterogeneous-column
references, exposure composition, downstream heating, very thick transmission,
whole-call budgets, late source mismatch and cancellation.

## Prescribed exterior states with inviscid entropy admission

`MixtureEuler1d::step_entropy_consistent_with_exterior_states` extends the
original first-order, fixed constant-area SSP-RK2 owner to prescribed
species-bearing exterior gas at either physical end. Exterior velocity may be
nonzero. These are fixed Riemann traces using the original variable-caloric
entropy-stable interface flux; the existing stationary compact-reservoir
characteristic closure retains its separate semantics. An absent exterior
requires a stationary reflecting end. Two absent exteriors execute the original
closed entropy entry exactly. Raw cells retain the caller's responsibility for
canonical source-species order and energy references.

For mathematical entropy `eta=-S`, each end supplies
`-v_b dot Q_b + sign*dt*A*psi_b`, where `Q_b` is the actual inward SSP-integrated
species/momentum/referenced-energy transfer and
`psi_b=u_b*sum_i(R*rho_bi/M_i)`. The effort includes momentum, kinetic energy and
the advective entropy potential. It is cached from ORIGINAL caloric interval
operations in the same immutable reference basis as the existing endpoint
entropy reducer; source formation or entropy datums are not subtracted as large
independent energy/entropy densities.

Before the native accepted-field/clock swap, the actual original-to-final state
change minus exterior supply must satisfy one absolute entropy allowance.
Each species separately balances its actual boundary transfer, scaled by
original inventory plus absolute throughput at both stages and both ends.
Only the explicitly bounded numerical continuity reference defect is removed.
The report retains source-state change, signed source entropy supplied at each
end, admitted generation and the reference defect. Rejected attempts contribute
no accepted transfer. Final cancellation and all refusal paths retain accepted
state/time.

Common positive species support, or common absence, is required at every face.
Unsupported reconstruction, moving/variable-area geometry and unowned open
ends refuse. This does not establish a general positivity/invariant-domain
theorem, finite-reservoir feedback, source-composition admission, exact moving
contact equilibrium or experimental flow validity. Focused regressions include
independent variable-Cp SSP state/transfer references, moving exterior entropy,
source-datum invariance, a reflecting end, common absence, closed-path identity
and atomic cancellation.

### Complete reacting entropy with prescribed exterior states

`ReactingEuler1d::step_entropy_consistent_with_exterior_states` composes the
original C(dt/2)-A(dt)-C(dt/2) owner with the prescribed-exterior entropy flow
entry. It retains the original reactor, chemistry method, accepted/trial
workspaces and bounded whole-composition restart loop. If A shortens the common
duration, both chemical halves restart from the original accepted state; only
the accepted middle transport's physical boundary inventory can enter the
complete account. Both absent exteriors delegate to the closed entropy entry.

The native SSP owner retains opaque actual per-end conservative transfers,
signed time-area measures and absolute per-species stage/end throughput.
The complete reacting account contracts those transfers again in its ORIGINAL
pre-chemistry entropy basis. It does not sum component entropy reports. With
`S=F+sum(alpha_i*m_i/M_i)`, actual chemical changes `Q_i` and exterior transfers
`B_i`, admission is
`Delta F + sum(alpha_i*Q_i/M_i) - I_boundary_normalized >= -tolerance`.
Every species independently balances `Delta m_i-Q_i-B_i`, with absolute
chemical and boundary throughput included in its continuity scale. Real
chemical source entropy is retained; boundary and reaction transfers are
removed at the species-residual level before applying the entropy reference.

All original source-state inversions, temperature restrictions, chemical
element/charge/energy checks, final conservative reporting and cancellation
precede the original accepted-field/clock swap. Refusal exposes no accepted
prefix. Cached exterior efforts and one reusable optional receipt slot allocate
only on open calls; retries reuse them. Existing closed molecular/RKL interval
paths retain their original arithmetic and scope.

This entry is one accepted inviscid reacting step with fixed prescribed
exteriors, first-order space and constant area. It does not enable molecular
half steps, finite exterior depletion, moving geometry, complete open reacting
intervals or a general positivity/order/combustion-validity theorem. Focused
regressions cover variable Cp and nonzero species entropy constants, independent
complete state-minus-boundary reconstruction, accepted-duration retry replay,
chemical-plus-boundary continuity, closed dispatch and late cancellation.

### Open reacting molecular transport with complete entropy admission

`ReactingEuler1d::step_entropy_consistent_with_dissipation_and_exterior_states`
and `step_entropy_consistent_with_stabilized_dissipation_and_exterior_states`
connect the ORIGINAL source-conjugate molecular, chemistry and inviscid owners
in D-C-A-C-D with the same prescribed stationary exterior states. Both molecular
halves retain their Maxwell--Stefan, heat-conduction, viscous and supplied
constitutive cross-effect fluxes and their native entropy/source gates. No
separate reaction or viscous heat source is added. Moving molecular exteriors
refuse under the existing stationary-bath model; the inviscid-only exterior
entry retains its separate moving-trace support.

Each completed molecular half exports only already-admitted actual signed end
species/momentum/energy transfers and double-double absolute stage throughput.
Its advective time-area measure is zero: molecular exchange does not receive an
Euler entropy-potential term. The whole composition contracts D1, A and D2 in
one immutable pre-D source entropy basis. It admits actual final state entropy
minus ALL external supply while retaining real chemical entropy. Individual
species satisfy reaction-plus-boundary continuity with both absolute chemical
and boundary throughput. Component entropy reports remain available but are
not summed into complete-state admission.

Original accepted/trial storage, common-duration restart and final temperature,
conservative-report and cancellation checks remain authoritative. Rejected or
shortened attempts clear all physical exchange receipts. RKL2 retains one
actual RHS counter across both halves and every discarded whole attempt;
work exhaustion remains terminal. Its signed boundary quadrature and absolute
negative-weight throughput remain owned by the original recurrence. Exporting
an accepted inventory introduces no fallible operation after molecular
publication and no additional EOS or constitutive evaluation.

Two absent exterior states delegate to the original closed SSP/RKL paths and
map their exact reports to zero-boundary form. This entry covers one accepted
first-order, fixed constant-area reacting step with prescribed stationary baths;
it does not establish finite-bath feedback, a full-duration open window,
unconditional positivity, global temporal order or combustion validation.
Regressions cover variable Cp/nonzero species entropy references, independently
reconstructed complete boundary supply, both-half and discarded-attempt work
limits, closed compatibility, moving-bath refusal and atomic cancellation.

## Frequency-resolved photoionization with evolving opacity

`chemistry::photoionization::spectrum::advance_spectral_photoionization`
evolves constant, bandwidth-integrated photon groups through one fixed-area,
uniformly mixed cell. Every band retains its original cross section and its
own incoming, absorbed and transmitted photon and energy accounts. Changing
neutral abundance changes each band's absorption; an average photon energy
or average opacity is not substituted for the spectrum.

Exact monochromatic subflows compose symmetrically. Adaptive full/two-half
comparison controls estimated second-order local errors in population,
electron temperature and individual band accounts. Cross sections prepare
once; trials reuse bounded buffers without allocating or solving a nonlinear
system. Positive subflows retain finite neutral tails and photon budgets.
Explicit source-temperature, group-count, step-count, charge, photon and
energy allowances gate publication, with cancellation during group work.

Each absorbed photon creates one ion and electron. The original chemical gap
is stored once, and the band's excess energy enters the changing electron
population's sensible energy. Initially absent electrons need no artificial
seed. The accepted endpoint retains the original thermodynamic bank for
electrical, collision and optical consumers without an LTE projection.
Seven focused regressions cover the exact monochromatic limit, independently
integrated two-band heating, temporal and band-order refinement, neutral
tails, empty/dark/subthreshold inputs, initially unseeded mixed spectra,
conserved reference shifts and transactional budget/cancellation refusal.

The scope is discrete beams through a fixed mixed cell, with quasi-static
propagation and immediate photoelectron thermalization. It does not infer
spectral quadrature accuracy, spatial multigroup fronts, secondary ionization,
recombination during the exposure, or empirical/global error bounds.

## Spatial spectral ionization and evolving shielding

`chemistry::photoionization::spectrum::column::advance_spectral_column`
propagates the original discrete photon groups through ordered, individually
mixed cells of a common beam area. Every cell retains its own atomic source,
thermodynamic bank and actual populations. Each complete monochromatic column
subflow passes the integrated escaping photons into the next cell using the
original depletion/transmission kernels. Its fluence solution accounts for
changing upstream opacity; downstream fluence-equivalent means are not treated
as simultaneous constant illumination during a multigroup exposure.

Symmetric composition of those whole-column band flows resolves concurrent
spectral hardening. Adaptive full/two-half comparison controls local estimated
population, electron-temperature and individual band-transfer errors. Prepared
cross sections and reusable band-major buffers keep trial work linear in the
cell/group product without per-trial allocations or nonlinear solves. Explicit
cell, group, cell/group-product and step limits bound work and storage.

Per-cell incoming/absorbed/outgoing band accounts close across each interface;
per-band and whole-column reductions independently account for exiting light,
chemical storage and electron sensible energy. Actual midpoint/end source
temperatures, original charge admission, energy allowances and cancellation
gate one complete publication. Inputs remain unchanged on every refusal.
Focused controls cover the monochromatic-column and single-cell limits,
independent coupled ODE evolution, downstream hardening, temporal refinement,
interface/global conservation, dark bands and transactional refusal.

The model retains fixed geometry, quasi-static propagation, local mixing and
immediate photoelectron thermalization. Spatially distinct cells resolve
shielding but do not resolve structure within a cell. Recombination during
transport, finite light-travel time, scattering, secondary ionization and
spectral-quadrature or empirical accuracy remain separate capabilities.

## Radiative capture with escaping photon energy

`chemistry::radiative_recombination::advance_radiative_recombination`
advances the original atom/ion/electron bank at fixed volume using reciprocal
photoionization-derived rate and captured-energy moments. Each event removes
one ion and one electron, creates one neutral, and emits one escaping photon.
The photon receives the original chemical energy plus the captured electron's
kinetic energy. Heavy temperature and spectator inventories remain unchanged.

The electron energy account includes the changing electron population.
Preferential capture can raise the temperature of the remaining electrons
while electron sensible energy and total gas energy decrease. Chemical
release is not deposited into the electron pool. Conserved elemental/charge
formation references cancel before evaluating the gas-plus-radiation balance.

A midpoint prediction evaluates both moments at the same stage temperature.
The ensuing exact binary-reactant update preserves finite reactants and has
a continuous equal-population limit. Full/two-half-step comparison controls
the second-order method for smooth coefficients. Stages allocate no species
vectors; final publication uses the existing frozen ionized-gas state.
Every changed source temperature requalifies its numerical moments and
conditional tails. Explicit work, positivity, temperature and energy limits
refuse without mutating the inputs. Empty reactants and zero duration retain
the exact original state.

Regression cases compare coupled population, temperature and escaping-energy
evolution with independent DOP853 references, including increasing electron
temperature during radiative cooling, refinement, source-reference invariance,
exact depletion limits, cancellation and bounded-work refusal. The model
assumes immediate electron thermalization and escaping recombination photons;
it supplies no trapped-radiation transport, stimulated capture, complete
multilevel recombination or closed-gas entropy inequality.

## Simultaneous local collisional and radiative ionization

`chemistry::collisional_radiative::advance_collisional_radiative` evolves the
same original ground-state atomic channel with electron-impact ionization,
its three-body reverse, local photoionization and spontaneous radiative
capture competing during one interval. The reciprocal photon source retains
the same atomic bank and original rate and captured-electron-energy moments.
Prescribed local photon fluxes are externally maintained; this entry does not
attenuate a column or consume a finite photon inventory.

An exponential-midpoint update integrates a frozen positive two-population
production/destruction law and its energy moments. Actual electron density
and temperature refresh its coefficients, and full/two-half comparison
controls second-order local estimates for populations, reaction throughput
and temperature. Original-source evaluation, source bounds, step budgets and
cancellation remain explicit; nonlinear trials allocate no species vectors.
Integrated event and energy accounts use cumulative ODE-state error scales;
normalizing by each individual time-step increment would impose an artificial
extra step restriction. This positive frozen-rate update does not make the
remaining temperature-dependent reaction feedback unconditionally stable.

Gross photon-driven creation and radiative capture remain distinct from net
chemical change. Absorbed photon energy enters once, escaping recombination
photons remove chemical plus captured kinetic energy, and electron sensible
energy follows the changing number of electrons. A neutral initial state can
acquire physical photoelectrons without an artificial seed. The final state
retains its original thermodynamics and separate heavy/electron temperatures
for material and electrical consumers, without an LTE projection.
Actual published population differences independently gate the event account,
and the gas-plus-radiation balance uses those actual endpoint differences.
Focused cases compare all four competing channels and both radiation accounts
with independently integrated closed-form moments, recover the analytic
fixed-temperature Riccati limit from zero electrons through repeated turnover,
and check refinement, formation-reference invariance and atomic refusals.

This is a rigid, uniformly mixed, single-ground-state channel with fixed heavy
temperature, immediate electron thermalization and escaping emitted photons.
It does not supply electron/heavy thermal exchange, excited-state kinetics,
stimulated emission, radiation trapping or a closed-gas entropy inequality.
The rate-equation basis is NIST's FLYCHK manual, section II.B.2;
<https://nlte.nist.gov/FLY/Doc/Manual_FLYCHK_Nov08.pdf>. Source-model uncertainty,
conditional rate-tail bounds and estimated temporal errors remain distinct.

## Conservative isothermal gravity with evolving mixture density

`multicomponent::IsothermalMixtureGravity` prepares an original-species
isothermal atmosphere on a stationary constant-area mesh.
`MixtureEuler1d::step_with_isothermal_gravity` advances the existing species,
momentum and energy coordinates with reflecting ends and first-order or MUSCL
spatial reconstruction. The reference retains its original species/caloric
bank, canonical mass fractions, physical mesh, temperature, pressure datum and
signed acceleration. Species densities remain dynamic; density is their
compensated sum, including when local composition changes.

The existing `try_new` constructor retains its uniform-composition mechanical
atmosphere. `try_new_diffusive` instead interprets the supplied fractions at
the pressure datum and uses every original molar mass:
`rho_k(x)=rho_k(x0)*exp(M_k*g*(x-x0)/(R*T))`. Heavy and light species therefore
have different scale heights; total pressure, energy and heat capacity come
from the original species EOS. At fixed temperature, each species' chemical
potential per mass plus gravitational potential is constant. This reference
assumes separately conserved neutral ideal-gas species; active charged species
are refused because no electrostatic balance is supplied.

Preparation caches original source-admitted reference face states and
tractions, and analytic exponential cell averages evaluated stably to floating
precision. The native Riemann and wall kernels remain the flux owners.
Subtracting the same reference tractions before differencing preserves
reference rest, while the density perturbation supplies the actual buoyancy
force. The report separates physical gravity impulse from the represented
reference-traction balance defect.

The diffusive profile integrates each partial density independently. A stable
expm1 expression retains the zero-gravity limit, with logarithmic product
evaluation only when intermediate exponential overflow or underflow would
hide a representable result. An authored zero species stays exactly absent;
an unrepresentable positive species is refused, never silently removed.

MUSCL reconstructs density and pressure relative to the integrated reference
cell averages, with physical velocity and the existing common mass-fraction
limiter. Both actual physical-face traces are admitted through the original
composition-dependent caloric source before accepting a slope. The same
contraction applies to every reconstruction coordinate, retaining linear
composition constraints and exact-zero species without clipping or reseeding.
An unchanged three-cell reference stencil has exactly zero perturbation slopes;
a neighboring disturbance still enters the original conservative Riemann solve.
The original first-order path remains the native fallback, and the higher-order
path reuses the existing slope and stage storage without additional allocation.

For the diffusive reference, composition is reconstructed as a departure from
the reference fractions. The physical and reference fraction gradients use
the same common limiter; their difference retains a zero fraction sum to
roundoff. Adding the physical reference-face fractions preserves the original
composition axis without renormalization. Both faces still pass the native
EOS admission and one common slope contraction. Large departures that would
produce an inadmissible fraction or seed an actually absent species use the
original positive-density, composition-preserving first-order fallback.

Gravitational work uses the same species mass fluxes as each native SSP-RK2
stage. The potential convention is phi=-g*(x-x0), with the explicit discrete
cell-center sum phi_i*M_i. Internal face transfers cancel in the combined gas
and gravitational energy account. This is a discrete potential quadrature,
not the exact within-cell integral of the continuous atmosphere. Per-cell
potential changes and source exchanges use compensated reductions; represented
source-energy loss remains a diagnostic, not an invented heat source.

Original source identity, formation-energy references, EOS/caloric admission,
CFL checks, retries and final cancellation remain authoritative. Energy and
potential accounts are evaluated before the physical state/time swap. Failed
stages remain private scratch. Zero gravity delegates the original unforced
entry, preserving its bits. Preparation owns the added buffers; successful
steps add no scratch allocation.

The scope excludes tapered/moving meshes, open reservoirs, self-gravity,
finite-rate species diffusion and a jointly advanced chemistry model.
Preparing a diffusive-equilibrium reference does not supply a transient
molecular-diffusion solver, nor does it relax arbitrary species inventories
toward that reference. The broad conservative gas control-volume bead remains
open.

The MUSCL regressions exercise exact stratified rest on a stretched mesh for
both reconstruction orders, real composition-driven buoyancy, source/formation
invariance, species and gas-plus-potential conservation, original CFL/caloric
retries, late cancellation and zero-gravity bit parity. Independent uniform
free fall checks actual acceleration, temperature and composition; an integrated
analytic small-pressure perturbation checks second-order momentum-response
convergence under mesh refinement. This is not a general convergence or entropy
theorem.

The diffusive-reference regressions additionally cover a 100:1 molar-mass
ratio, independently integrated partial densities and source energy,
chemical/gravitational potential balance, exact rest for both spatial orders,
small-pressure perturbations at two amplitudes under mesh refinement, actual
species/energy conservation, formation-gauge invariance, reference-present
species absent from the actual flow, representational refusal and final
cancellation. The standalone production profile leaf has independent
antiderivative, partition-additivity, small-gravity and overflow/underflow
checks; it uses the same deterministic elementary-function owner as the flow.

Four new integration tests cover analytical reference rest, composition-driven
buoyancy, species and combined-energy conservation, formation-gauge invariance,
native caloric/CFL retries, late cancellation, source/geometry admission and
zero-gravity parity. The new production file and corrected integration tests
passed native rustfmt. Root and an independent reviewer checked the discrete
force/work and atomic-publication scheme.

Focused native validation is complete for the exact published source in
`57475f27b2af915c77f30f9c6eb8a3166289258e`. After rebuilding genuine current
thermochemistry, the original chemistry/body-force/dynamics/multicomponent
source hierarchy compiled successfully. All four new gravity tests passed.
Three unchanged native regressions also passed: accepted-work integration,
zero-force bit parity, and invalid-source/late-cancellation atomicity. Every
compiler and test process returned exit code zero; no product correction was
needed after the interrupted run resumed.

This bounded native source check excludes unrelated fusion and optional
material root modules. It is not a whole-workspace Cargo or DSR pass. The
focused integration target remains
`cargo test -p fs-gas --features experimental-euler --test mixture_gravity`.

### Complete open reacting durations in one original entropy basis

`ReactingEuler1d::advance_stabilized_entropy_consistent_with_exterior_states_for`
and its `_to` companion complete the requested duration or represented target
through the ORIGINAL RKL2 reacting interval owner. Every accepted D-C-A-C-D
step uses the same stationary prescribed exterior states. Two absent states
delegate to the original closed entropy interval and map its exact reports.

The window retains its original thermochemical source, caloric domain,
reference temperature, species reference densities and entropy constants.
Local steps collect fresh initial temperatures and inventories while borrowing
that immutable basis and cached exterior efforts. No per-step report rebasing,
source-entropy reevaluation or sum of entropy allowances is used to complete
the window. Source/domain, reference-basis and exterior-effort compatibility
are required before merging a complete admitted step.

Only accepted chemical changes and actual molecular/advective boundary
transfers enter the compensated window account. Chemistry, boundary entropy,
species transfers and absolute stage throughputs merge transactionally.
The final gate reduces the original conserved rollback rows and actual final
state directly, retaining real chemical entropy and applying ONE whole-window
allowance. Per-step generation and species residual extrema are diagnostics.

The native interval keeps one reused conserved backup, shared coupled-attempt
and molecular-RHS work limits, actual final-remainder integration and exact
represented clock semantics. Any late entropy, source, work or cancellation
refusal restores the whole original field and clock, including an accepted
internal prefix. A zero duration still admits the original operator/source,
stationary baths and entropy scope.

The fixed-area, first-order, stationary prescribed-bath and common species
support restrictions remain. These entries do not evolve finite exterior
storage, establish higher-order positivity, or supply a global temporal-error
estimate. Six focused regressions cover independent complete source/boundary
entropy reconstruction, work exhaustion after an accepted prefix, late
cancellation and reuse, zero-duration admission, represented clock refusal,
and rejection when individually tolerated local entropy decreases exceed the
single original-window allowance.

## Reacting compression of the common material parcel

With `material`, `material::compression::advance_closed_parcel_compression`
advances the original `ClosedGasParcel` through a finite-time prescribed volume
stroke while its original mechanism reacts. Volume follows a constant
logarithmic rate from the current to the requested volume. Each caller-bounded
interval applies an original-material isentropic half-volume map, original
finite-rate rigid chemistry, then another isentropic half-volume map. Product
composition, NASA calorics, pressure, transport and acoustic properties remain
owned by the existing material and chemistry sources. Formation heat is already
in those calorics and is never added a second time.

This symmetric composition has second-order splitting accuracy for exact
subflows; native chemistry error remains controlled by its selected integrator.
The caller must refine `substeps` to assess compression/reaction coupling error;
there is no automatic splitting-error certificate. Work is the sum of original
isentropic caloric increments, with positive work delivered to the gas on
compression. The report retains stored energy change, energy-minus-work defect,
source-stage entropy change and chemical step counts. The energy defect is
compensated directly from the chemical material-reconciliation residuals, so
large work totals cannot erase a small reported defect by subtraction. The
entropy diagnostic sums original material volume-step changes and native
reactor changes before material endpoint reconciliation; it excludes that
reconciliation's entropy correction and does not independently bound the
published parcel's total entropy change. This is thermodynamic volume work,
not an independent pressure quadrature or piston force measurement.

Requests are bounded to 1..=4096 trajectory intervals. Fixed-volume evolution
uses one original reactor call; zero duration cannot move volume. Invalid,
nonrepresentable, cancelled, out-of-domain or failed intermediate states never
publish partial parcels. No constant-gamma surrogate, immediate equilibrium,
piston dynamics, shock, wall heat transfer or phase-change model is implied.
Manufactured regressions cover adiabatic compression/expansion laws, an analytic
exothermic reacting-compression solution and second-order refinement, independent
endpoint energy/work closure, original rigid-reactor compatibility, late
cancellation and source-domain refusal. They are numerical/physical coupling
checks, not empirical kinetic validation.


### Adaptive reacting material compression

`material::compression::adaptive::advance_closed_parcel_compression_adaptive`
adds caller-owned absolute
temperature [K], original species amount [mol] and interval work [J] tolerances,
plus a relative allowance on actual interval changes. Each trial evaluates one
full interval and two half intervals on the original prescribed logarithmic
volume trajectory. The SAME private original-owner interval composition serves
fixed and adaptive calls. Original species order is bound once; allocation-free
ordered comparisons include initially absent products without trace seeding.
Carrier inventories, arbitrary formation-energy origins and work accumulated
before the trial do not enlarge the local error scales.

Acceptance uses raw full/fine discrepancies in all three physical channels,
without an assumed-order divisor or Richardson extrapolation. Only the
two-half-step material state and its physical work, energy, residual and stage
entropy enter the accepted trajectory. The native chemistry integrator retains
its own accuracy, source, conservation, entropy and per-call work controls.
These local comparisons do not certify a global solution error or complete
published-parcel entropy change. Symmetric exact-subflow composition follows
Strang, SIAM Journal on Numerical Analysis 5(3), 506–517 (1968),
<https://doi.org/10.1137/0705041>; it does not make approximate native subflows exact.

Explicit minimum/initial/maximum time steps, at most 4096 accepted comparison
intervals and at most 12288 coarse/fine interval evaluations bound adaptation.
All rejected and coarse trials count against the evaluation budget; reports
separate their chemical work from accepted fine physical totals. Distinct
represented full/half/quarter times are required. The exact final remainder is
integrated and the exact requested terminal volume/time retained. A minimum-step
accuracy failure, exhausted budget, native physical/source refusal or late
cancellation preserves the entire original parcel, including after an internally
accepted prefix. Native physical refusals remain terminal; only comparison
failures trigger smaller trials. Fixed volume and zero duration retain the
original single-reactor route. External heat remains zero.

Focused manufactured regressions cover independent reacting temperature/work
solutions during compression and expansion, work-driven refinement, collision
chemistry whose composition error is independent of caloric changes, accepted
and rejected work bounds, exact endpoints, source-domain and late cancellation
rollback, and fixed-volume compatibility. The exponential analytic fixture
explicitly declares irreversible first-order A-to-B chemistry; the shared
reversible material fixture remains unchanged.

## Reacting electron and heavy thermal pools with radiation

`chemistry::collisional_radiative::advance_collisional_radiative_with_thermal`
extends the existing simultaneous impact/three-body/photoionization/radiative
capture owner with original neutral and Coulomb heat exchange, independent held
electron/heavy powers and evolving heavy temperature. The basic entry retains
its fixed-heavy, zero-exchange behavior. `CollisionalRadiativeThermal` carries
the original collision sources, signed powers, heavy temperature support and
caloric inversion bound; the result retains those inputs.

The same `PreparedExchange` owner used by existing gas energy evolution queries
actual trial populations, volume and both temperatures. No collision coefficient
is reused at an unrelated state. The original atomic channel has equal constant
monatomic neutral/ion molar heat capacities; its heavy nucleus count is
conserved, while spectator amounts remain fixed. Original-source heavy sensible
energy differences and capacities therefore support the shared bounded caloric
inversion without double-counting the chemical gap or imposing a constant
mixture heat capacity on spectators.

The midpoint stage freezes conductance and capacities and integrates the
linearized thermal difference exponentially. Its forcing includes changing
electron population as well as chemical/radiative energy and external work.
Constant capacities, conductance, populations and powers recover the exact
linear exchange solution; general temperature-dependent reaction/transport
feedback retains the local midpoint estimate and is not unconditionally stiff
stable. The full/two-half comparison includes both temperatures and the
accumulated transfer/work accounts. Trials allocate no species vector.

Electron-to-heavy transfer appears once with opposite signs in the two pool
balances. Neutral and Coulomb contributions remain distinct. Published pool
and total balances independently use actual endpoint populations and original
source calorics, including final population roundoff. Pure thermal/power
evolution proceeds even when chemical rates vanish. A nonzero held electron
power at exact initial electron absence refuses; photoionization with zero
held electron power can create the first carriers and activate exchange.
Source/domain exits, exhausted work or caloric inversion, and cancellation
publish no partial endpoint. The model remains a rigid ground-state Maxwellian
parcel with escaping emitted photons and immediate electron thermalization.

Focused regressions compare neutral thermal relaxation in both directions with
its independent nonlinear closed form, and the full reacting/radiating parcel
with independently integrated analytic Maxwellian/Milne/Landau moment equations
(DOP853 and Radau references). They also exercise variable spectator calorics,
formation-reference shifts, pure power with inactive chemistry, initial neutral
photoelectron creation, original-source admission and cancellation. Numerical
agreement on these manufactured laws is not experimental atom/plasma validation.


## Photon transport with competing reactions and thermal evolution

`chemistry::radiation_transport::advance_radiation_transport` couples ordered
mixed gas cells to a finite incident photon exposure. Each cell retains one
original reciprocal atomic source for photoionization, electron impact,
three-body recombination and spontaneous radiative capture. Recombination can
restore neutral opacity during illumination; each later photoionization
consumes another incident photon. Gross photon/capture traffic remains separate
from net pair creation. Optional per-cell `CollisionalRadiativeThermal` inputs
use the same local thermal kernel and heavy caloric inversion described above.

Symmetric composition interleaves local reaction/thermal midpoint subflows with
exact monochromatic whole-column attenuation/depletion. The spectral bands use
their own symmetric composition. A common adaptive full/two-half comparison
controls the combined local second-order estimate in populations, both thermal
pools, reaction integrals and band photon accounts. Original cross sections and
canonical source bindings prepare once. Reused contiguous point and band
buffers avoid trial state cloning, trial species allocation and nested adaptive
public solves. Accepted plus rejected attempts and cell/group extents are
explicitly bounded; no global, spatial or empirical accuracy is certified.

Per-band incident/absorbed/transmitted counts close at every cell interface and
the whole-column boundary. Photon totals retain double-double accumulation and
direct transmission evaluation instead of subtracting nearly equal incident
and absorbed counts. Published actual gas energy closes against separately
accumulated beam absorption, escaping recombination photons and both external
pool works. Internal electron/heavy exchange cancels in the column account.
The absolute energy allowance applies to the column as well as each cell; it
is not silently multiplied by cell count. Empty columns transmit the original
exposure, darkness still permits chemical and thermal evolution, and exact
electron absence requires no artificial seed.

Geometry and beam area are fixed and propagation is quasi-static. The emitted
recombination radiation has explicitly unit escape probability: its energy
leaves the system and it is not fed back as diffuse ionizing radiation. This
does not model photon trapping, light-travel delay, scattering, stimulated
emission, excited states, hydrodynamic motion or a non-Maxwellian distribution.
The supported source interval must contain the numerical split stages; a
temperature-domain trial failure can reduce the step, while other original
source and work refusals remain explicit. Failure leaves all caller cells and
sources unchanged and publishes no partial column.

Five focused regressions cover an independent simultaneous two-cell reaction
and attenuation ODE; repeated photoionization and capture beyond the original
neutral inventory; second-order refinement in a resolved lower-density fixture
with unchanged optical depths; combined beam, neutral exchange and both pool
powers; and darkness, empty/zero-time identities, bounded work and cancellation.
The dense adaptive reference remains separate from the resolved fixed-step
refinement fixture. Exact conservation alone does not establish temporal order
or validate the declared unit-escape approximation.

## Primitive transport coefficients and separated face-energy contributions

`PlanarTransport::coefficients_from_primitive` and
`SpatialTransport::coefficients_from_primitive` query the original transport
provider at an admitted cell primitive. They use the original caloric trace,
temperature, pressure and mole fractions, reuse binary-coefficient scratch, and
validate the returned shear/bulk/thermal and binary transport coefficients.
They do not construct a face, solve a Maxwell--Stefan system or advance time.
Cancellation and provider errors remain explicit.

`TransportFaceReport` also carries `bulk_viscosity_pa_s`,
`nonmechanical_energy_flux_w_m2` and `nonmechanical_diffusivity_m2_s`.
The separated energy contribution retains the original species-enthalpy,
Fourier and reciprocal Dufour terms before mechanical stress work is added.
The corresponding diffusivity excludes the mechanical stress estimate.
Engineering, entropy-conjugate and wall/bath paths populate these receipts;
the original face flux keeps its existing arithmetic. This permits a caller
to replace the molecular stress discretization without subtracting a huge
old stress work from a small physical heat or enthalpy flux. The separated
receipts do not add an independent energy source or a new entropy guarantee.

## Reacting material compression with finite-rate heat exchange

`material::compression::thermal::MaterialThermalCompressionReactor` adds held-bath
and finite-source-body heat exchange to the same common `ClosedGasParcel` volume
stroke. It prepares one original `MaterialIsochoricReactor` and reuses its chemical
stages and dense feed workspace throughout and between strokes. The fixed path
scheduler and symmetric isentropic-volume / rigid-reacting-thermal / isentropic-
volume composition are shared with the original adiabatic compression owner;
the fixed and adaptive adiabatic APIs retain their existing behavior.

The middle stage integrates reaction and signed `G (Tpartner - Tgas)` together
inside the original implicit chemistry/caloric solve. NASA temperature-dependent
Cv, formation energy, composition, pressure, source transport and acoustic
response remain original. There is no frozen gamma, nominal heat capacity,
separate temperature update or second reaction-heat term. The caller supplies
constant total conductance in W/K; no moving-area or material-name conductance
law is inferred.

A held bath has fixed temperature and infinite capacity, with finite declared
conductance. A finite `FiniteHeatBody` instead carries its exact original
enthalpy/phase source, mass, current enthalpy, density and held body pressure
through all accepted intervals. Original Cp variation and latent plateaus remain
inside the existing finite-body stage elimination. Body volume changes never
replace the prescribed gas volume. Source and neutral dilute-gas admission,
caloric bracketing, thermal/chemical entropy gates, mass conservation and final
material reconciliation remain with the original owners.

The gas report separates pressure work delivered to the gas, signed partner heat,
actual stored-energy increments and `Delta Ugas - W_on_gas - Qgas`. Pressure work
is the original isentropic caloric work, not an independent measured-pressure
quadrature. Gas first-law residuals are accumulated directly from the material
rigid stages, so subtracting large stroke totals cannot hide their small defects.
The finite-body result additionally retains its original `Delta H`, `Delta U`
and outward pressure-volume work. Its two directly accumulated combined residuals
represent `Delta Ugas + Delta Hbody - W_on_gas` and
`Delta Ugas + Delta Ubody + W_by_body - W_on_gas`. Exchanged heat cancels once;
body expansion work and gas compression work retain distinct signs and owners.

Reported gas entropy sums original volume-map and native chemical/thermal-stage
changes. Partner entropy is `-Q/Tbath` or the finite body's original integral of
`dh/T`; reported production sums native gas-plus-partner production and the
volume-map increments. As with adiabatic compression, native chemical entries
precede material endpoint reconciliation and exclude its entropy correction.
These stage diagnostics do not independently bound published-parcel entropy or
certify global passivity, integration accuracy or experimental validity.

`MaterialThermalCompressionSolve` admits 1..=4096 equal-time intervals on the
logarithmic volume path. The exact requested final volume and elapsed duration
are retained. Fixed volume uses one original reactor call regardless of the
requested interval count; zero duration preserves both physical states after
admission and cannot move volume. All proposed gas/body states remain private
until the entire stroke and final cancellation check succeed. Intermediate source,
transport, balance, representability, work-budget or cancellation refusal rolls
back the complete stroke, including an internally accepted prefix. A failed
proposal changes reusable numerical workspace only. Native integration controls
and their energy/entropy tolerances apply separately to each stage; accumulated
residuals are reported, not relabelled as a global error bound. The caller must
refine splitting intervals to assess coupled heat/reaction/compression error.

Focused manufactured regressions cover a closed-form exothermic bath trajectory
through both compression/cooling and expansion/heating, second-order splitting
refinement of temperature and pressure work, independently integrated variable-Cv
and finite variable-Cp body trajectories, latent melting and both work accounts,
original fixed-volume and insulated compatibility, exact no-op states, late
cancellation/replay and source refusal after a feasible prefix. These fixtures
are numerical physics checks, not measured kinetic or heat-transfer calibration.
Pinned Rust parsing/formatting was checked; the new native tests have not been
executed because the shared host has no available compilation memory or disk.


## Coupled atomic charge-state networks with shared electrons

`chemistry::ionization_network::advance_ionization_network` advances multiple
original atomic charge-stage ladders through one shared electron population and
temperature. Each original edge supplies impact ionization, its thermodynamic
three-body reverse, photoionization and reciprocal radiative capture. Sources
must share the exact original thermodynamic bank and structural electron.
Canonical lower/upper species ordering makes input permutations equivalent.
Duplicate reaction identities or edges, branching stage definitions and cycles
refuse; disconnected ladders and omitted frozen spectator stages are explicit.

One symmetric composition couples positive exponential-midpoint edge updates.
Every edge queries the current shared electron density and temperature. A
single full/two-half comparison controls species, electron temperature and each
edge's integrated reaction and radiation accounts. Prepared indices and reused
population/account buffers avoid per-trial species allocation and nested public
adaptive solves. Accepted plus rejected attempts, channels, and all edge rate
queries have explicit bounds. This is a local second-order split estimate, not
an unsplit, unconditionally stiff, global-error or empirical accuracy claim.

Publication independently reconstructs each actual species change against
integrated stoichiometry, each connected nucleus inventory, and shared charge.
Original adjacent energy gaps define levels relative to each conserved ladder
root, removing conserved elemental/charge reference offsets before subtraction.
Actual endpoint level storage and electron sensible energy then close against
separately accumulated absorbed and escaping photon energy. Heavy translational
population roundoff is included at the fixed heavy temperature. Per-edge gross
photoionization and radiative capture remain distinct from net collisional and
net charge-stage progress. Exact neutral darkness stays neutral; illumination
may create the first electrons without a seed. No partial state publishes after
cancellation, source refusal, exhausted work or a failed conservation check.

This lane holds geometry, heavy temperature and local incident intensities
fixed. The maintained field supplies an external photon reservoir; its finite
transport and diffuse recombination feedback are not supplied by this owner.
Emitted photons escape and photoelectron kinetic energy thermalizes immediately.
Excited states, molecules, charge exchange, heavy-pool exchange/evolution,
continuum lowering, hydrodynamics and electrical work are not inferred.
Focused tests compare a simultaneous independent analytic-source H/He ladder
ODE, resolved second-order refinement, canonical source ordering, conserved
formation-reference shifts, first photoelectron birth, and transactional
admission/work/cancellation cases. NIST level gaps with manufactured cross
sections establish numerical checks, not validated H/He plasma kinetics.


## Diffuse spectral recombination radiation in mixed slabs

`chemistry::diffuse_recombination::evaluate_diffuse_recombination` solves the
frozen, quasi-static bound-free radiation field of ordered homogeneous
plane-parallel slabs with common area and vacuum boundaries. Original Milne
spectra supply isotropic spontaneous emission; original energy-dependent cross
sections supply absorption. Exact uniform-volume-source characteristics include
both self-absorption and bidirectional nonlocal illumination. The source spectrum
is never replaced with one mean photon. No beam, scattering, excited-level
cascade, stimulated emission, light-travel delay or radiation-pressure work is
inferred.

Positive angular Gauss quadrature uses mu=t^4 to resolve grazing rays. Spectral
panels split at every original source knot, onset, thermal scale and individual
emission cutoff. Separate Gauss 8/16 estimates check every cell and directional
interface number/energy moment. Original total Milne moments independently
check the integrated emission. Each emitter's omitted tail retains a conditional
original-envelope bound, propagated conservatively to absorption and interface
observables. Cold emitters may omit a bounded high-energy emission tail while
their original high-energy opacity remains active. Shared finite absorption
support, unresolved tails, exhausted angular/spectral/total work or cancellation
refuse a complete field rather than extrapolate or publish partial results.

The result retains N+1 bidirectional interfaces, source states, gross capture
and self/nonlocal absorption, escaped photon energy, net ionization and
electron/chemical/gas power rates. Actual photon and independently transported
energy balances close in every cell and globally; source/sink internal transfers
cancel. Trial quadrature reuses prepared source/geometry bindings and contiguous
scratch buffers without cloning states or banks. This query does not evolve
populations or temperatures; its nonzero rates are inputs to a later coupled
evolution owner. Quadrature differences are estimates, conditional tail bounds
are not physical-source uncertainty certificates, and no spatial or experimental
accuracy is implied.

Focused verification covers independent isotropic slab escape from optical
depth 1e-6 through 1e6, two-slab nonlocal angular integrals and reflection, an
independently integrated cold/hot spectral-opacity case, uniform-slab subdivision
invariance, and original support/work/cancellation/darkness behavior.

## Prescribed reacting compression and expansion histories

`fs_gas::material::compression::history` adds
`advance_closed_parcel_compression_history` for a sealed reacting parcel driven by
explicit elapsed-time/physical-volume knots. The first knot binds time zero and
the exact current parcel volume. All volumes must be finite and positive, times
must increase strictly, and adjacent volume ratios must be representable.
Admission permits 2–4097 knots and checks every knot before evolution.

Volume is log-linear on each positive-duration segment. Equal adjacent volumes
are physical-time dwells through the original rigid reactor. Turns and dwells
remain explicit; returning to the initial volume does **not** collapse a history
to fixed-volume chemistry. Every segment retains the preceding accepted
candidate's actual species inventory, material query axes, mass and original
thermochemical/transport sources. There is no species reseeding, source
reconstruction from a catalog approximation, or second kinetics solver.

The typed solve mode delegates to the existing owners:

- Fixed mode uses the declared original substep count on each moving segment and
  one original reactor interval per dwell. Its separate whole-history substep
  allowance is in 1–4096 and is never replenished at a knot.
- Adaptive mode retains the original temperature, original-species inventory and
  interval-work comparisons. Its maximum accepted intervals and maximum native
  evaluations apply to the **entire history**, including coarse and rejected
  probes. Each dwell costs one interval/evaluation. A known dwell may consume a
  final allowance of one or two evaluations despite the original adaptive
  control admission minimum of three; moving comparisons still require all
  three evaluations before starting. The preceding native step suggestion seeds
  the next segment within the unchanged native step bounds and survives dwells.
  No comparison crosses a prescribed rate discontinuity.

The report chronologically compensates the original accepted work, internal
energy change, independently accumulated reconciliation residual and stage
entropy. It separately counts all native probes/chemical work and retained
physical intervals, so coarse and rejected trajectories do not enter physical
balances. The reported end time is the exact final knot; a compensated sum of
segment-reported durations minus that knot exposes floating-point subtraction
rounding between absolute times. This diagnostic is not a time-integration error
certificate. Native adaptive accuracy/work/time errors carry the global history
clock and consumed prefix counts; source failures retain their enclosing segment
bounds without inventing an unreported partial native clock.

The complete caller parcel is published only after every segment, aggregate
balance validation and final cancellation checkpoint succeed. Invalid later
knots, exhausted global allowances, native accuracy/domain failures and late
cancellation leave the complete original parcel and source receipts intact.

Focused manufactured regressions cover an exothermic compress–dwell–expand
cycle against the independent solution
`dT/dt = a*T + Q/(3R)*4*exp(-4t)`, with second-order fixed-step refinement and
net cycle work. An equal-caloric `2A -> A+B` mechanism independently checks
`N_A(t) = 1 / (1/N_A(0) + 4*integral(dt/V))`: species remain sensitive to
pressure history even when final volume and temperature return to their initial
values. Other checks exercise exact global budgets including the last
single-evaluation dwell, translated failure clocks, independent binary-duration
roundoff, source identity, and rollback after an accepted prefix.

This is a prescribed homogeneous sealed-parcel history, not a piston-motion or
compressible-flow solver. It adds no heat exchange, rate calibration,
extrapolation or empirical combustion claim. The original source-domain,
chemical-conservation, material-energy, stage-entropy and local adaptive error
limitations remain in force; accumulated local comparisons are not a global
accuracy or total-entropy certificate.

The new Rust sources were parsed/formatted with the pinned toolchain. Native
compilation and these regression tests remain unexecuted because the shared
workspace has no available compilation memory or disk; no passing result is
claimed here.


### Direct net radiation rates in opaque slabs

The diffuse transfer owner integrates net reaction, electron-energy and gas-energy
rates directly from nonlocal absorption minus locally escaping emission.
Self-reabsorption cancels algebraically before rounding. Subtracting separately
integrated gross emission and absorption can otherwise erase a finite escape
loss in a strongly opaque slab. Gross source, absorption and directional
interface accounts remain independently available and checked. Signed net
moments participate in the same spectral, angular and conditional tail controls.
The net rates also close against actual incoming minus outgoing interface
throughput, with roundoff scaled to that traffic rather than gross self-recycling.

The uniform-slab regression extends through optical depth 1e20 and compares net
recombination, gas cooling and captured-electron cooling directly with the
independent finite-escape result. This establishes numerical behavior of the
declared frozen slab model; it does not add a dense-plasma or time-dependent
radiation approximation.

### Real-gas nozzle momentum and stationary side-wall traction

The original NASA/Peng--Robinson field now admits stationary varying-area
passages through `RealGasEuler1d::new_on_mesh`. Both reconstruction choices and
the existing force/wall compositions use the same physical `Mesh1d` volumes.
The geometric pressure source is `p_i (A_right - A_left)`; it exchanges axial
momentum with the stationary passage and performs no energy work. The shared
ideal/real-gas update pairs this source with each face pressure before area
multiplication, so a uniform-pressure resting state is preserved without
subtracting large extensive pressure forces. Finite-speed steady nozzle states
are not exactly well balanced.

Both SSP-RK2 stages accumulate the original cell-pressure side-wall force with
compensated summation. `StepReport::area_momentum_impulse` integrates those
stage forces, and the inlet/outlet changes use their own physical face areas.
Mass and energy remain boundary/source balances; area pressure enters only
momentum, once. The initial step estimate, face CFL gates and both mechanical
force endpoints use both adjacent area/volume metrics. Source forcing and
finite-wall heat/drag retain their separate impulse/work/heat accounts.

No geometry re-integration or extra per-step allocation is introduced.
Varying-area volumes remain the geometry producer's physical integrals;
constant-area volumes retain the axial consistency check. A periodic seam
requires identical end areas at construction and boundary replacement, with
transactional refusal. The same original material branch, energy gauge,
source admission, retry and cancellation limits apply. This is the smooth
quasi-one-dimensional engineering model, excluding abrupt-junction mixing,
moving cross sections, general viscous stress and a fully discrete entropy
theorem.

Five native tests of the actual shared geometry kernel cover pressure rest,
moving-flow advection, extensive mass/momentum/energy balance, second-order
smooth pressure-source quadrature and scale-robust face CFL arithmetic.
Five original-material integration regressions are included for stationary
balance, closed compression, forced area scaling, periodic seam admission and
cancel/retry publication. Those full real-gas integration tests were not run
in this constrained environment; the focused native result does not replace
the repository Cargo/DSR gate.

### Adaptive heat exchange during reacting material compression

`material::compression::thermal::adaptive` adds
`advance_with_heat_bath_adaptive` and `advance_with_finite_body_adaptive` to the
prepared thermal compression reactor. Every full/two-half comparison uses the
same original `advance_interval_with` composition, native simultaneous reacting
thermal stages and original volume-work owner as fixed thermal compression.
The prepared kinetics workspace is reused; each probe owns its gas state and
finite-body enthalpy/phase state. The two fine intervals carry the first half's
actual body state into the second. Coarse and rejected states publish nothing.

Five independent raw comparison channels control gas temperature, every original
species amount, gas pressure work, exchanged heat and actual finite-body enthalpy.
The last channel compares original body mass times the full/fine specific
enthalpy difference, including latent inventory when temperature is unchanged.
Specific enthalpies are subtracted in compensated precision before mass scaling.
Relative allowances use interval changes, so formation/enthalpy datums, large
carrier inventories and earlier work/heat cannot dilute the error criterion.
The original ordered species comparison is shared with adiabatic adaptation;
canonical species keys are bound once, including initially absent products,
without cloning source cards or rebuilding per-trial maps.

Only accepted fine work, heat, stored-energy changes, body expansion work and
original source residuals enter physical totals. Chemical work across all
successful probes, including coarse and rejected trajectories, is separately
counted with checked arithmetic. Gas and coupled gas/body first-law residuals
remain direct compensated sums of their original-owner residuals; they are not
reconstructed by subtracting large aggregate energies. Body heat is never added
again to gas formation heat or pressure work. Existing source-stage entropy
diagnostics retain their documented material-reconciliation limitation.

Positive absolute allowances, a relative allowance in [0,1), ordered finite time
bounds, at most 4096 accepted comparisons and at most 12288 original interval
evaluations are explicit. Every moving full/fine comparison costs three probes.
Rejected work consumes the same global evaluation budget. Native chemistry
retains its own per-call work, accuracy, mass, energy and entropy gates.
Full, midpoint and both quarter times must be distinct represented values; the
exact final remainder and requested terminal volume are retained. Fixed volume
uses one native reacting thermal interval after adaptive-control admission;
zero duration preserves the original gas and body exactly. Whole-stroke
accuracy, work, source/domain, representability or final-cancellation refusal
preserves both original endpoints, including after an accepted private prefix.

Step selection reuses the original integrator-dependent heuristic. Acceptance
always checks the raw discrepancy with no Richardson divisor, extrapolation,
species clipping or nominal heat capacity. These local comparisons include
native integration differences; tolerances below that accuracy may refuse at
the minimum step. No global time-error bound, unconditional stiff accuracy,
piston inertia, spatial wall law or physical validation is claimed.

Focused regressions cover analytic bath compression/cooling and expansion/
heating, heat-only refinement, body-enthalpy-only refinement through a latent
plateau with shifted energy datums, density-dependent collision inventory
control, original variable gas Cv and body Cp, distinct body expansion work,
one-call fixed-volume compatibility, exact no-op states, accepted/rejected work
limits, minimum-step and represented-clock refusal, late cancellation/replay
and source refusal after a feasible prefix. Native execution and pinned
formatting have not run because the execution environment is unavailable;
source/API review and independent analytic reference calculations are separate
from runtime verification.

## Thermal reacting compression histories and finite-body cycles

`material::compression::history::thermal` extends the prepared
`MaterialThermalCompressionReactor` with
`advance_history_with_heat_bath` and
`advance_history_with_finite_body`. Both use the same admitted
`CompressionKnot` schedule as the original adiabatic history entry:
2–4097 knots, exact initial volume binding, positive finite volumes, strictly
increasing elapsed times and representable adjacent volume ratios. Each segment
has its own logarithmic volume path. Turning points and physical-time dwells are
retained even when the last volume equals the first.

The existing prepared thermal owner supplies every fixed or adaptive endpoint.
The history layer adds no chemistry, volume map, heat-transfer stage, caloric
inversion or time controller. Original gas species amounts, formation energy,
source query axes, neutral-gas admission and transport/caloric sources continue
through every knot. The held bath retains its declared temperature and total
conductance. A finite body retains its original source curve, mass, held
pressure, conductance and current enthalpy/phase, including heat-flow reversal
on a latent plateau where its temperature is unchanged. Its density-derived
volume never substitutes for prescribed gas volume.

Fixed subdivisions apply separately to moving segments, while one explicit
whole-history allowance bounds their total; each dwell costs one native
interval. Adaptive mode delegates to the original thermal full/two-half
controller, preserving its temperature, species, work, heat and finite-body
enthalpy comparisons. Accepted-interval and evaluation caps apply to the entire
history, including coarse/rejected probes. The last one or two available
evaluations can admit a known one-evaluation dwell; a moving comparison still
requires three. A native suggested step carries across knots and dwells, within
the same original bounds. No comparison crosses a prescribed rate discontinuity.

The generic history report retains the original gas or gas/body physical
payload. Chronological compensated sums preserve work delivered to the gas,
signed partner heat, gas storage, independently accumulated first-law residual,
gas/partner stage entropy and production. Finite-body histories additionally sum
the original body enthalpy/internal-energy changes, outward pressure-volume
work and both coupled residuals. Exchanged heat cancels once; neither chemical
formation energy nor latent storage becomes an extra heat source. All-probe
chemical work remains separate from accepted physical accounts. An explicit
duration discrepancy reports compensated segment durations minus the exact
final elapsed knot.

Native adaptive refusals retain the global history clock and spent prefix
counts in the typed thermal-adaptive history error. Complete knot admission
precedes evolution. Gas and body remain private proposals through every
segment, aggregate-account check and final cancellation checkpoint; a late
work limit, source-domain exit, accuracy failure or cancellation publishes
neither endpoint. Numerical reactor workspace may change on refusal and remains
reusable for deterministic replay.

Focused G1/G3/G4 regressions cover a closed reacting bath cycle against an
independent exact first law with fixed-step refinement; a variable-Cv gas and
variable-Cp/density body against independent coupled balances; latent melting
and partial reversal through both fixed and adaptive histories; exact global
budgets including a final single-evaluation dwell; translated error clocks;
source refusal after a feasible prefix; late cancellation/replay; and an
independent binary-duration rounding check. Manufactured irreversible
equal-caloric cases stop before their mixing-entropy maximum rather than
loosening native entropy gates.

This remains a prescribed homogeneous gas and lumped thermal partner. It adds
no spatial wall gradients, piston inertia, moving-area conductance inference,
gas phase equilibrium or measured combustion/heat-transfer calibration. The
original native source-domain, material-reconciliation entropy and local-error
limitations remain unchanged. Aggregating local comparisons and residuals does
not certify global integration accuracy or coupled passivity.


## Direct radiative energy accounting for evolving ionization

The shared collisional-radiative reaction kernel retains absorbed-minus-emitted
photon energy as a directly accumulated signed moment. Its published
`net_radiation_energy_j` and electron/total energy checks use this moment rather
than subtracting gross cumulative emission and absorption. Gross photon events,
absorbed energy and captured-electron energy remain separate physical accounts.
The finite incident-beam owner retains the same signed moment in its cell and
whole-column reports while keeping its independent finite-beam balances.

A crate-private collision/thermal subflow excludes both photon processes when
a spatial diffuse owner supplies the common radiation field. Signed diffuse
increments change the original adjacent populations and electron sensible
energy together, preserve the original heavy temperature and clock, and admit
first-electron birth from actual photoelectron energy without a trace seed.
Direct net event/energy closure is checked separately from gross consistency
at the gross accounts' floating-point resolution. An inadmissible trial refuses
instead of clipping populations or temperatures. The spectral transfer owner
must supply its independently conserved signed moments; this seam alone is
neither a radiation solve nor a new physical closure.

Focused G0 regressions retain a finite net loss below the representable spacing
of gross self-recycling accounts, verify first-carrier energy and transactional
refusal, and prove that the collision subflow excludes both photon processes.
The original simultaneous reaction/thermal and finite-beam regressions remain
the compatibility checks for the shared kernel.


## Thermal feedback in shared-electron charge-stage networks

`chemistry::ionization_network::advance_ionization_network_with_thermal`
extends the original canonical multielement charge-stage owner with evolving
electron and heavy temperatures, original neutral/Coulomb energy exchange,
and separately held electron/heavy powers. It consumes the existing
`CollisionalRadiativeThermal` inputs and prepared thermal/caloric kernel; it
does not fit a replacement collision law or call nested adaptive solvers.

Symmetric thermal half-steps surround each complete canonical reaction sweep.
Consequently heat exchange and external work act once per parcel interval,
independent of the number of active charge stages. Every thermal stage refreshes
the complete current population vector, including charge-squared screening and
energy exchange from multiply charged ions. Every reactive edge uses the
current shared electron temperature and heavy temperature. One full/two-half
controller includes both temperatures, all species and reaction accounts, and
the signed neutral/Coulomb transfer integrals. Original source quadrature,
caloric iteration and complete step budgets remain explicit.

All populated or reachable neutral targets must have an original collision
source before evolution begins when exchange is selected. The shared Coulomb
law retains its original bank/volume and weak-coupling/classical validity
checks at changing populations and temperatures. Heavy calorics integrate the
original models, including variable-Cv spectators; equal admitted monatomic
heat capacities across each reactive ladder preserve this original caloric
reference. Actual endpoint heavy-count roundoff is included in the energy
check. Independent electron/chemical and heavy balances close against their
own held work and equal/opposite exchange, and the total closes against
absorbed/emitted radiation plus both external works.

The fixed-heavy entry remains available. The thermal entry requires at least
one admitted atomic source to bind the electron caloric convention. A nonzero
held electron power refuses an initially absent electron pool; with zero such
power, photoionization can create the first carriers without a seed. Sources,
controls and original state are retained, and failure publishes no partial
endpoint. Volume and maintained local illumination are fixed, emitted photons
escape, and no diffuse multichannel radiation, transport, electrical-field
solve, molecules, excited levels or dense-plasma correction is inferred.

Focused G1/G3/G4 tests compare a simultaneous independent H/He/He++ chemical
and two-pool energy ODE with neutral and finite-mass Coulomb exchange, verify
second-order refinement and conserved formation-reference shifts, and exercise
neutral birth, missing-source/work/domain refusals and cancellation. Cross
sections are manufactured test sources; these checks do not validate an
empirical multielement plasma model or certify global temporal error.


## Prescribed ionizing beams with diffuse recombination radiation

`chemistry::diffuse_recombination::evaluate_diffuse_recombination_with_illumination`
admits explicit discrete photon energies, directional fluxes and signed direction
cosines at either plane-parallel boundary. A ray's face throughput is its flux
times face area times the absolute direction cosine; its path optical depth is
the normal column depth divided by that cosine. Original cross sections are
resolved once into a bounded source table. Exact homogeneous-slab attenuation
uses `expm1` for absorption so optically thin losses survive when transmission
rounds to one. Subthreshold photons pass without ionization or heating.

The incident field adds to the existing spectral/angular diffuse solution.
Incident absorption, self-recycling and photons received from other slabs remain
separate accounts. All three feed the same net ionization, original chemical
storage and photoelectron excess-energy moments. Every cell and the full column
close emission plus incidence against absorption plus outgoing photons/power.
An empty column transmits both directions unchanged. Positive beams must lie
inside every original absorber's support, including initially absent absorbers;
invalid directions, source support, work/allocation or cancellation refuse
without changing caller state. The vacuum entry point is unchanged.

G1/G3/G4 tests compare independent Beer--Lambert attenuation for two slabs,
opposite oblique beams and different photon energies, exact frozen-field
superposition, subthreshold transparency, absorption below transmission
roundoff, empty columns and transactional refusal. These are numerical tests,
not validation of the supplied atomic data or a spectral/angular discretization
error bound for the caller's illumination. The formal transfer equation and
beam/diffuse decomposition are consistent with
<https://journals.ametsoc.org/view/journals/atsc/57/10/1520-0469_2000_057_1635_rtiwim_2.0.co_2.xml>;
ionization and heating by transported recombination photons are discussed in
Raičević et al. (2014), <https://arxiv.org/abs/1311.0182>.

## Evolving diffuse recombination radiation and gas feedback

`chemistry::diffuse_recombination::evolution::advance_diffuse_recombination`
evolves the actual populations and electron energies of ordered mixed slabs
while transporting their original spontaneous recombination spectra. Each slab
retains one original adjacent atomic channel, fixed volume and common face
area. Optional original neutral/Coulomb collision sources and held powers evolve
the heavy pool through `CollisionalRadiativeThermal`. Emission changes the gas,
changed populations alter opacity, and the evolving electron temperature
changes subsequent emission and collision rates. Boundary escape, self
reabsorption and illumination of other slabs are distinct physical accounts.

The frozen query and live owner share the same spectral/angular transfer kernel
and its conservation gates. A prepared workspace retains source/geometry
bindings and reuses cell, observable, spectral-grid and angular buffers across
raw population/temperature evaluations. No trial clones a thermodynamic bank
or frozen state, creates a species vector, or nests a public adaptive solve.
Bounded field evaluations and total radiation work include predictor and
rejected-trial work; each field retains its spectral/angular/tail controls.

Symmetric local collision/thermal half-steps surround an explicit-midpoint
diffuse update. Both photon processes are excluded from those local half-steps
to avoid double-counting emission or absorption. One full/two-half temporal
controller includes populations, both temperatures, reaction and energy
integrals, and every directional interface exposure. Positivity and original
temperature support are acceptance conditions: an inadmissible radiation
trial reduces the step instead of clipping or creating a trace carrier floor.
This is a jointly estimated second-order composition, not an L-stable method
or an unconditional positivity guarantee for arbitrary steps.

The transfer kernel supplies direct signed net ionization, electron heating
and gas-energy moments. Local self-recycling cancels before rounding, retaining
escape loss even when independently rounded gross emission and absorption are
equal. First-electron birth uses actual absorbed excess energy. Integrated
photon and power fluxes cross the same slab interfaces as the reaction/heat
updates. Publication independently closes actual original gas storage plus
outer-boundary radiation loss against external pool work, and net radiative
pair creation against escaped photons. The absolute gas-energy tolerance
applies to each cell and once to the complete column. Caller states and sources
remain unchanged on source, budget, conservation or cancellation refusal.

This is quasi-static transport in transversely infinite plane-parallel slabs
with vacuum boundaries and immediate electron thermalization. It adds no
incident beam, finite light-travel delay or stored photon inventory, gas
motion, scattering, stimulated emission, excited-level cascade, radiation
pressure or multiple competing absorbing channels within a slab. The separate
shared-electron network retains its maintained local field and unit-escape
closure. Numerical time and spectral estimates do not certify spatial
resolution, model-form error or experimental plasma accuracy. For physical
context on diffuse ionization and heating, see Raicevic et al., MNRAS 437,
2816 (2014), https://arxiv.org/abs/1311.0182; that different transport method
is not a validation oracle for this implementation.

Focused G1/G3/G4 cases compare an independent analytic-Milne/E3 evolving-slab
ODE, nonlocal first-electron birth independent of an empty pool's stored
temperature, and simultaneous impact/reverse capture, diffuse transfer,
neutral thermal exchange and held work. Resolved timestep refinement checks
second-order behavior. A formal optical-depth-1e20 long-exposure case retains
finite net recombination and electron/gas cooling while gross self-recycling
exceeds 1e17 mol. Empty/dark/zero-duration and original-source, allocation/work
admission, step and cancellation cases retain transactional behavior. These
manufactured numerical tests do not validate those extreme parameters as a
physical plasma or bound the quasi-static approximation.


## Competing atomic absorbers in one gas mixture

`chemistry::ionization_network::transport::evaluate_network_radiation_transport`
transports prescribed discrete boundary beams through ordered, uniformly mixed
cells containing multiple original atomic charge-state ladders. Every cell uses
one shared species vector and electron pool. At each energy the normal channel
optical depths add; exact uniform-slab attenuation consumes one shared photon
supply, and each channel receives its opacity fraction of the absorbed photons.
This is the species-partition rule in Pawlik and Schaye, arXiv:1008.1071,
Section 4.1. Co-located species are not ordered as separate shielding layers.
Positive/negative direction cosines select the entry boundary; projected face
flux and oblique optical path are both represented. The small absorbed fraction
uses `expm1`, independently of transmitted flux rounding.

Original source identity, bank, charge ladder and supported photon energies are
admitted for all channels, including initially absent but reachable stages.
Cross sections and canonical channel indices prepare once, and radiation stages
reuse fixed buffers. Explicit cell/channel/beam products bound allocation;
source preparation, characteristic sweeps and radiation evaluations have work
limits and cancellation checks. Result rates retain the original state and
canonical sources, per-channel ionization and excess-energy deposition, and
per-beam local and whole-column photon accounts. No population or temperature
change is implied by the frozen query. Empty and transparent cells transmit.

The internal network radiation adapter applies simultaneous edge increments
through their original stoichiometry and a single shared-electron caloric update.
It rejects combined stage depletion without clipping, and retains a compensated
signed radiation-energy source independently of gross emission/absorption.
Actual chemical/sensible storage still has the original absolute energy audit;
the separate gross-traffic comparison allows only its floating-point roundoff.
Geometry is fixed, propagation is quasi-static, and the caller owns discrete
spectral/angular resolution. This entry does not supply excited-state channels,
scattering, secondary nonthermal electrons, diffuse emission or light-travel time.

## Diffuse radiation shared by co-located charge-stage networks

`chemistry::ionization_network::transport::diffuse::evaluate_network_diffuse_recombination`
retains each original mixed cell and its canonical reciprocal atomic sources.
Every emitting edge supplies its own Milne spectrum; all bound-free absorbers
in that slab compete through the sum of their actual optical depths. Exact
uniform-source characteristics carry isotropic emission across a common-area
plane-parallel column. Same-source recycling, absorption by another channel
in the same cell, nonlocal absorption and prescribed boundary-beam absorption
are distinct accounts. The existing finite-beam kernel is reused.

Cross-channel recycling transfers the difference between the original chemical
thresholds to the shared electron pool. Pairwise local exchanges are assembled
antisymmetrically before angular transport, while independent compensated cell
moments retain net reaction, electron energy and gas energy before large
per-channel exchanges are rounded. Trapped escape is therefore not obtained
by subtracting gross emitted and absorbed traffic. Cell and boundary photon,
power and chemical/electron balances are checked separately.

Every active emitter and absorber must support the common spectral interval.
Cold emission cutoffs do not switch off higher-energy absorption. Source
conditional tails and spectral/angular estimates are bounded independently of
incident beam brightness. Unsupported opacity, inadequate finite source tails,
invalid original-state/source bindings, cancellation and exhausted work or
quadrature budgets refuse the whole result. Retained source tables and scratch
buffers avoid per-channel state clones inside spectral/angular quadrature.

This entry evaluates a frozen field. It does not assume excited-state cascades,
case-B branching, a fitted escape fraction, an average emitted photon,
scattering, photon travel time or atomic source data beyond the originals.

## Optically-thin free-free cooling in atomic charge-state evolution

`advance_ionization_network_with_optically_thin_bremsstrahlung` and
`NetworkRadiationCell::with_optically_thin_bremsstrahlung` attach the original
`ElectronIonBremsstrahlungV1` to the shared atomic electron pool. The exact
source bank and structural electron must match. At least one admitted atomic
edge establishes the original Maxwellian electron caloric convention. Constant
or original Itoh nonrelativistic integrated Gaunt laws retain their source
temperature, charge, nondegeneracy, relativity and stationary-ion limits.

Each cooling stage refreshes the actual electron and charge-squared ion
densities and the temperature-dependent Gaunt factor. In the coordinate
sqrt(Te), a source-evaluated midpoint freezes C=P/sqrt(Te) for the interval.
The step is exact for a constant Gaunt law at frozen populations and second
order otherwise. Stable direct energy loss avoids subtracting nearly equal
thermal stores. Predictor and endpoint admission reject source-boundary or
zero-temperature crossings. Cooling surrounds the complete shared reaction
sweep, not each channel, and therefore applies once independently of edge count.

The native adaptive comparison includes escaped continuum energy and the
direct combined radiation source. Initial source admission and every rejected
query consume the shared rate budget. Bound-free photon counts and
`net_radiation_energy_j` retain their original meaning.
`free_free_escaped_energy_j` records only continuum escape;
`net_gas_radiation_energy_j` independently retains bound-free deposition
minus that escape before separately published components are rounded. Strict
electron and gas storage audits use this direct combined source. Finite-beam,
diffuse and thermal owners share the same cooling implementation.

This is an explicit optically-thin escape closure for the free-free component,
including when bound-free photons are separately transported or recycled.
Neither a spectrum, free-free photon count, absorption/reprocessing nor an
escape factor is inferred from a bolometric Gaunt source. The bare-Coulomb-ion,
weak-coupling and negligible collective-effect scope of that source remains
the caller's physical model.

## Evolving mixed ionization networks with diffuse recycling

`advance_network_diffuse_recombination` shares the native network time
controller with finite-beam evolution. Each radiation midpoint uses the
current original populations and electron temperature to refresh reciprocal
emission, common opacity, cross-channel recycling and nonlocal transport.
The surrounding native collisional/thermal steps include impact ionization,
three-body capture, elastic pool exchange and retained external powers;
spontaneous capture occurs exactly once through the diffuse field.

All channel increments update the common species vector together. Independent
compensated cell totals determine electron number, electron energy and signed
gas radiation energy after consistency checks against the canonical edges.
This preserves the small net loss of a trapped mixed cell without losing
large opposite chemical exchanges. Exactly neutral cells can acquire
electrons from transported photons; their birth temperature follows the
actual deposited photoelectron energy, without a seed population or a
stored-empty-pool temperature prescription.

The existing full/half-step comparison controls populations, temperatures,
channel ledgers, beam exposures, same-source/cross-source/nonlocal diffuse
traffic, direct signed cell moments and every directional interface
exposure. Rejected trials contribute to one global work/evaluation budget.
All source and scratch storage is retained across stages. Invalid trials
refuse or reduce the step through the existing controller, never clip
populations or publish a partial endpoint.

`NetworkRadiationResult::diffuse_radiation` exposes accepted interface
exposures and independent photon/source/boundary audits; the cell accessor
exposes spontaneous and recycling accounts. Existing beam-only entry points
and band balances retain their semantics. Gross emitted radiation is
distinguished from actual boundary escape. The strict actual-gas versus
direct-source energy check retains the original absolute chemistry
tolerance; only incoming/outgoing boundary subtraction receives its explicit
floating-point traffic allowance.

The field remains quasi-steady and plane-parallel. Time-step estimates and
spectral/angular/tail bounds are separate controls, not a global accuracy
certificate or an unconditionally stiff solver. Excited-state cascades,
scattering, moving radiation boundaries, finite light-travel time and
unsupported atomic cross sections are not inferred.

## Evolving mixed ionization networks under finite beams

`chemistry::ionization_network::transport::advance_network_radiation_transport`
uses that same prepared shared-opacity field at every radiation stage. An
explicit-midpoint update applies all source increments simultaneously to one
population vector and Maxwellian electron pool. Thus a newly created He+ stage
can absorb in a later stage of the same exposure, and changing upstream mixtures
alter downstream illumination. Exact absence of electrons is supported using
the combined physical photoelectron energy, without a seed or channel-order
temperature artifact. Both boundary directions and oblique incidence remain
supported throughout the evolution.

Symmetric composition surrounds the common photon step with the existing dark
network's impact ionization, three-body capture, escaping radiative capture,
and optional original thermal exchange/held work. No maintained local photon
field is also supplied to that kernel. One full/two-half comparison controls
populations, both temperatures, original channel ledgers, direct signed energy
and every integrated incoming/absorbed/outgoing beam count. All trial buffers
and source tables are reused; global radiation work and rate-query budgets
include rejected trials. Negative populations and unsupported trial temperatures
reduce the step; they are never clipped. No unconditional stiffness, global
time error or empirical accuracy claim follows from the local estimate.

Publication independently checks channel absorption against transported photons,
preserves source-retained endpoints and canonical channel accounts, and audits
actual chemical/electron/heavy storage against direct radiation and held work.
Gross radiation-energy comparison has an explicitly reported roundoff allowance;
this does not relax actual gas conservation. Empty columns transmit their whole
exposure, zero duration preserves identity, and failures publish no partial
result. Volumes and geometry are fixed and propagation is quasi-static. All
spontaneous recombination photons escape this entry; it does not yet transport
diffuse recombination between multiple co-located absorbing charge stages.

## External illumination during diffuse ionization evolution

`chemistry::diffuse_recombination::evolution::advance_diffuse_recombination_with_illumination`
evolves the original atomic populations and electron/heavy thermal pools while
the same prescribed discrete boundary beams pass through the changing opacity.
Every explicit-midpoint radiation stage combines beam attenuation with the
current spontaneous spectrum and its self/nonlocal reabsorption. The existing
collision/thermal half-steps still exclude both photon processes; neither
incident work nor diffuse emission is counted twice. Source cross sections are
prepared once and retained across all accepted and rejected trials.

Incident photon absorption and its energy have distinct per-cell integrals.
Their net reaction and excess-energy effects enter the original chemical and
electron accounts. Integrated incoming and outgoing interface traffic closes
against actual gas storage and held external pool work. A neutral illuminated
column cannot take the dark shortcut: the first carriers are created from
absorbed photons with their actual excess energy, not an arbitrary seed or the
temperature stored on an absent electron pool. Empty columns preserve their
complete transmitted exposure; zero duration preserves identity.

The final gas-storage versus direct-radiation-source audit retains the absolute
chemistry energy tolerance. Boundary-energy subtraction alone additionally
allows `128*f64::EPSILON*(incoming_energy+outgoing_energy)`, reported as
`boundary_energy_roundoff_j`, because the directional interface rates are f64:
a nearly transparent bright beam may have equal rounded incoming/outgoing
energies while `expm1` still resolves physical absorption. The separately
reported `gas_source_residual_j` cannot use this throughput allowance to hide
a gas-heating error. The allowance is an arithmetic admission tolerance, not
a certified transport-error bound. A live thin-depth regression independently
checks retained ionization/heating and this separation of the two audits.

The existing temporal controller now includes incident absorption integrals
alongside all population, temperature, diffuse and interface accounts. Source
support, cancellation, field work, step and conservation failures publish no
partial result. The fixed-volume, mixed-slab, quasi-static assumptions remain:
beams are held over the requested interval, no photon flight time or stored
radiation field is introduced, and the caller supplies the discrete spectral
and angular approximation. This does not resolve excited channels, scattering,
gas motion or competing absorbing charge stages inside one slab.

New G1/G3/G4 regressions use the independent self-shielding solution
`exp(a*n_neutral(t))-1 = (exp(a*n_neutral(0))-1)*exp(-a*b*t)` in a formally
dilute column where original recombination is negligible; compare first-carrier
temperature against absorbed excess energy; retain simultaneous incident and
diffuse recycling with independent energy/charge checks; and exercise empty
exposure, zero duration, work and cancellation refusals. These new tests were
authored after the execution environment became unavailable and have not run.
No native build, temporal-order validation or empirical accuracy is claimed for
this extension in the absence of executed results.

### Original-source real-gas moving chambers and mechanical forcing

With `material,experimental-euler`, `RealGasEuler1d::step_with_moving_walls`
advances a resolved, constant-area real-gas field between two prescribed
impermeable ends. Both speeds use the increasing axial coordinate and remain
held over the accepted interval. Affine face motion reuses the original
nonuniform `Mesh1d` geometric-conservation law. Extensive SSP-RK2 updates
mass, laboratory momentum and laboratory total energy with the same face
velocities and stage volumes. The first Euler stage, twice-displaced second
Euler endpoint and final average all resolve through the original
NASA/Peng--Robinson source before publication.

`RealGas::moving_hllc_flux` constructs relative-frame traces from original
density and temperature, retaining source-owned internal energy instead of
subtracting large laboratory kinetic energies. The admitted HLLC/HLL flux is
transformed back to laboratory conserved transport F-wU. Signal speeds remain
relative to the face, and both stage geometries receive their actual CFL gates.
`moving_reflecting_wall_flux` reuses the original shock/isentrope pressure
relation at fluid velocity relative to the wall. Numerical mirrored ends
retain their original source-checked Riemann pressure. Both kinds have exactly
zero mass transfer and pressure-times-wall-speed energy transfer. The same
stage fluxes determine the reported left/right impulse and pressure work.

`step_with_moving_walls_and_body_force` surrounds this transport with the
existing held-force half steps. The first uses accepted density and physical
volume; the second uses transported density and displaced physical volume.
Each endpoint retains original calorics, Mach admission and full-interval
face-relative CFL. Boundary pressure work and external source work remain
distinct and are counted once. Retry reconstructs all candidates from accepted
fields and geometry; cells, primitive caches, mesh and clock swap only after
the final cancellation check. Zero motion reuses the existing stationary
method exactly; zero force reuses the source-free moving method.

The numerical scratch and geometry buffers are allocated lazily and reused.
MUSCL uses the existing cached original-source traces; boundary cells retain
their original relation. This is a bounded second-order engineering lane on
smooth admitted solutions, subject to source limiting and first-order retries.
It does not infer piston inertia, finite wall heat, moving area, open advective
ends, a phase flash or an entropy/positivity theorem. The prescribed-wall lane
does not compose the separate finite-wall thermal/drag owner.

Seven executable native tests of the production kinematics module pass:
independent laboratory F-wU transport, pressure work, Galilean composition,
bitwise uniform translation, geometric conservation, nonuniform extensive
balance and second-order work refinement against an independent adiabatic
volume ODE. A metadata-only compile checks the changed production owner,
geometry, source and reconstruction control flow with explicitly unimplemented
constitutive facades; it is not an EOS validation or full crate build.
Seven original-material regressions in `tests/material_nonideal_flow/moving.rs`
specify physical moving compression, Galilean wall/face fluxes, uniform
translation, force-density work on displaced volumes, original-path identity,
formation-datum invariance and late cancellation/retry. Those complete
material integration tests were authored but not executed in this resource-
limited session; the focused arithmetic checks do not replace them.

### Pressure-driven real-gas piston inertia

With `material,experimental-euler`, `RealGasEuler1d::step_with_piston` couples a
fixed left wall and a finite-mass right piston through the original
NASA/Peng--Robinson ALE stages. Each midpoint-speed trial uses the pressure
impulse and pressure work from the selected numerical or physical wall law;
independent impulse and combined gas-plus-piston energy checks retain external
force work and exterior-pressure work. `step_with_suspended_piston` adds a
stationary linear spring and viscous damper, using actual midpoint displacement
and velocity while retaining spring storage and nonnegative damping loss
without injecting that loss as gas heat. Cells, cached thermodynamics, mesh,
clock and piston velocity publish together after complete source and mechanical
admission. This requires constant area, two impermeable ends, frozen composition
and the declared homogeneous stable real-gas branch. The bounded engineering
coupling does not compose wall heat, chemistry, distributed forcing, phase
change or additional moving masses, and does not establish a general FSI
stability or entropy theorem.

### Wall heat and finite storage coupled to piston work

`RealGasEuler1d::step_with_piston_and_walls` and
`step_with_suspended_piston_and_walls` surround the original pressure/inertia
solve with two wall heat/drag half intervals. The first source changes gas
pressure before the piston solve; the second uses actual transported density
and displaced cell volume. Independent local finite walls, shared material
wall inventories and their authored conduction links carry their first-half
storage into the second half. Both shared-source halves spend one attempt's
bounded source-work allowance. Original nonlinear or latent h/u wall laws and
their pressure-work convention remain intact; reservoir heat, prescribed
sliding-drive work, internal slip heating, piston pressure work, spring storage
and suspension loss retain separate accounts. Damping is not implicitly
returned as heat. Complete gas/piston/wall balances and post-heat full-interval
relative-face CFL are checked before gas, caches, geometry, clock, piston, walls
and reports publish together. Refusal restarts the entire interval from accepted
inventories; accepted wall bindings follow the moved mesh. Shared-link stage
receipts describe the second heat half, and the shared aggregate covers both.

Reusable source and geometry buffers avoid field clones inside midpoint-speed
trials. Backward-Euler heat keeps this symmetric source/pressure/source
composition first order. Conductances, wall amounts and sliding velocities
are explicit inputs; no changing contact area or thermal piston mass is
inferred. The existing constant-area, impermeable-end, frozen-composition and
homogeneous real-gas branch restrictions remain. Chemistry, distributed body
forcing and phase flash are not composed by these entry points.

### Distributed body-force work coupled to piston motion

`RealGasEuler1d::step_with_piston_and_body_force` composes authored acceleration
or force density with the original finite-piston pressure/inertia solve.
`step_with_piston_sources` additionally admits the existing suspension and wall
operator in the same transaction. Ordering is body half, optional wall half,
native pressure/inertia evolution, wall half, body half. The first source
endpoint is prepared once per attempted duration and reused by every
midpoint-speed proposal; the second kick uses transported density and actual
moved volume. Original mechanical kicks preserve internal energy and record
signed extensive impulse/work in `gas.source_change`; authored acceleration
does not infer a gravitational potential inventory.

Pressure-only mechanical admission remains independent. Complete gas, piston,
wall and spring balances include body work, reservoir heat, sliding drive and
damping exactly once. Both force endpoints retain Mach admission, and the final
state must satisfy full-interval relative-face CFL. Source refusal shortens the
whole interval without demoting MUSCL. Rejection/cancellation publishes no gas,
cache, mesh, clock, piston or wall state. Zero forcing follows the existing
entry without source allocation. Existing constant-area, impermeable-end,
frozen-composition and original caloric/phase restrictions remain;
backward-Euler wall heat keeps first-order composition accuracy.

Focused tests add independent acceleration/force-density impulse and kinetic
work oracles, same-step pressure feedback, finite wall/suspension balances,
late source refusal/retry and zero-force parity. These added tests have not
been executed in the current resource-limited environment.

### Two original real-gas chambers sharing a finite piston

`RealGasEuler1d::step_with_shared_piston` evolves two independent
NASA/Peng--Robinson wave fields through one finite inertial piston;
`step_with_shared_suspended_piston` adds the original stationary spring/damper.
Each midpoint-speed proposal calls the native moving-volume RK2 operator with
`[0,w]` on the left and `[w,0]` on the right. Both actual interface pressure
impulses and their conjugate gas works feed the existing scalar mechanical
solve. Both gas energies plus piston kinetic change balance the authored
external-force displacement work. The suspended entry also retains spring
storage and nonnegative damping, without heating either gas. Fixed outer
supports carry their original reactions and zero work.

The two fields retain independent source models, frozen compositions, energy
references and grids. Interface coordinate, constant area and accepted clock
must match exactly. A common acoustic/contraction limit covers both chambers
and the whole speed bracket. Refusal restarts both native fields; gas states,
cached primitives, meshes, clocks and piston publish atomically after the final
cancellation poll. Candidate buffers are reused without field/source cloning
per proposal. Thermal walls, body forcing, chemistry, phase flash and leakage
remain outside this adiabatic shared-piston entry; original source and near-wall
accuracy limits remain.

Five focused test groups cover unequal-grid equilibrium, bidirectional pressure
work and initial Newton acceleration, bilateral native moving-wall replay,
independent source-datum invariance with suspension, and late refusal/retry.
These added native tests have not been executed in the current resource-limited
environment.

## Physical-time thermal forcing during reacting compression

`material::compression::thermal::schedule` adds immutable `HeatBathSchedule`
and `ConductanceSchedule` programs to the prepared original material reactor.
Bath temperature and total conductance vary linearly between prescribed physical
time knots; a finite body's program changes conductance only. The first time is
zero, subsequent times advance strictly, and the last knot supplies exact stroke
duration. One zero-time knot is a no-op program. Admission bounds each borrowed
immutable profile to 4097 knots, validates its physical values and polls
cancellation. There is no callback-dependent forcing, extrapolation or inferred
moving-area, convection or spatial wall law.

Fixed and adaptive methods reuse the same native volume/reacting-thermal interval
kernel, original calorics and physical accounts. Each interval freezes forcing
at its own physical midpoint. Every forcing knot is an exact interval boundary,
so full/two-half probes cannot miss an interior heating/cooling excursion. A
varying program at fixed gas volume still enters the time comparison; an exactly
constant program delegates to the original fixed or adaptive endpoint API.
Fixed forcing knots split the requested uniform grid; the complete union is
admitted against the original 4096-interval bound before chemical evolution.
Adaptive accepted/evaluation caps remain global, including coarse and rejected
probes. Exact forcing remainders may be shorter than the declared minimum step,
but a failed comparison there cannot reduce below that minimum. Represented
midpoint/quarter clocks must advance, and the terminal time and volume are exact.

Finite-body programs bind the exact original initial conductance. The new
`FiniteHeatBody::with_conductance` operation changes that boundary coefficient
without reconstructing or changing original source, enthalpy/phase state, mass
or held pressure. Each accepted body carries the exact endpoint conductance;
midpoint values never escape as the terminal coefficient. Existing adaptive
body-enthalpy error control remains active through latent plateaus and shifted
energy datums. Gas heat, body enthalpy/internal energy, body outward pressure
work and gas compression work retain their original distinct accounts. Native
residuals remain direct compensated sums, with no extra reaction or latent heat.

All physical proposals remain private through the complete program and final
cancellation checkpoint. Budget, clock, source-domain, balance or cancellation
refusal preserves both caller endpoints, including after accepted internal
intervals or rejected probes. Bath entropy sums the native `-Q/Tbath` at each
interval's own sampled temperature; it is not total heat divided by the final
temperature. This approximates the varying-bath entropy integral and does not
model the external energy required to prescribe the bath itself. The existing
source-stage material-reconciliation and local-vs-global accuracy limitations
remain unchanged. Temporal refinement, original chemical tolerances and the
physical validity of the prescribed thermal program remain caller obligations.

Focused regressions use independent coupled balances for simultaneous bath and
conductance ramps during compression, expansion and rigid heating; fixed-step
refinement; heat-only adaptive control at fixed volume; original variable gas
Cv and finite-body latent inventory under changing conductance; shifted enthalpy
datums and exact terminal coefficients; constant/no-op compatibility; exact
forcing-knot and global work limits; unrepresentable clocks; late cancellation,
source refusal and deterministic replay.

## Original real-gas pressure-controlled discharge

With `material,experimental-euler`, `RealGas::pressure_outlet_flux` and
`RealGasBoundary::PressureOutlet` supply downstream-pressure boundary physics
on the retained NASA/Peng–Robinson source. A compression uses the source-owned
`material::nonideal_shock::pressure_compression`: at each density ratio,
`j² = (p2-p1)/(v1-v2)` and `e2-e1 = (p1+p2)(v1-v2)/2` close mass, momentum and
energy. The density Newton tangent follows the same EOS derivatives; the
existing bounded shock solver performs the full entropy, sonic and conservation
admission. No outer piston-speed solve or fitted gamma is introduced.

Expansion follows `s = constant` and integrates `du = -a d(rho)/rho` using the
same bounded logarithmic-density quadrature as physical wall waves. The pressure
root uses the original `dp/dln(rho)|s = rho*a²`. A fan spanning the face is
sampled at `u=a`; reducing backpressure beyond that point does not increase
discharge. The bracket grows locally and stops when pressure or sonic discharge
is reached, rather than requiring an arbitrarily cold remote tail. Fully
outgoing waves retain the interior flux, while sufficiently strong backpressure
can drive a compression shock into an initially supersonic outlet. Exactly zero
external pressure uses this same positive-density sonic limit; native EOS
evaluation never requires a zero-pressure or zero-density state.

`RealGasWallControl` supplies the explicit homogeneous source paths and work
bounds; its compression momentum tolerance also gates the outlet pressure
residual. A weak compression below the shock resolution uses the native acoustic
tangent, checks the actual isentropic source path and admits its pressure error
against that tolerance. Its front uses the mean native incoming-family speeds,
retaining first-order nonlinear shock motion at sonic and slightly supersonic
discharge. Pressure-only backflow is refused because the entering
entropy/thermal state is unspecified. Positive source pressure, density and
temperature bounds are retained throughout; no clipping, phase flash or global
classical-branch certificate is implied. This outlet is stationary; the existing
moving-domain entry continues to require impermeable end walls.

Fluxes transport the owner's common-reference energy and pressure work exactly
once, and enter both SSP-RK2 stages and the existing extensive boundary balance.
Late cancellation or source refusal publishes no partial field. Regressions
include the independent 75-digit dense PR shock anchor, dilute ideal expansion
and sonic limits, large formation-datum shifts, reversal, weak acoustic response,
two-stage field conservation and mid-step cancellation. The general shock
equations follow NASA TN D-8274, *Flow Through Normal Shocks*,
<https://ntrs.nasa.gov/citations/19770008402>; its nitrogen data do not calibrate
this authored binary PR model.

## Body acceleration in finite-rate reacting transport

`ReactingEuler1d::step_with_body_force` composes the original conservative
mixture source with finite-rate chemistry on the same species bank and energy
references. The accepted map is C(h/2), F(h/2), A(h), F(h/2), C(h/2).
Uniform or cellwise acceleration is multiplied by the actual mixture density
at each F stage; force per volume keeps its prescribed value. The existing
mechanical kick couples represented momentum increments to midpoint kinetic
work, so forcing supplies neither an independent heat source nor a second
chemical-energy datum.

`step_with_dissipation_and_body_force` adds the original molecular D(h/2)
steps on both sides. Viscosity, species diffusion and conductivity retain the
authored source-bound transport law. The original chemistry and molecular
integrators retain their own order; symmetric composition alone does not
establish a second-order certificate.

The original acoustic and whole-composition retries choose one accepted h.
If accelerated or heated transport shortens h, both chemistry halves, both
force kicks and both optional molecular halves restart from the accepted
state. Final caloric and authored Mach admission precede publication.
`transport.source_change` retains the actual mechanical impulse, work and
representation residual. Complete momentum and physical/gauged-energy
balances subtract its impulse and work once; the source representation
residual remains visible rather than being removed from conservation error.
Cancellation or refusal publishes no accepted cell or time change. Zero
forcing follows the original unforced numerical map.

The focused G0/G1/G3/G4 cases in
`multicomponent::reactive::tests::forcing` compare uniform accelerated heat
release with the independent homogeneous reactor and exact Newton impulse/
work; refine an independently known accelerating reaction/diffusion Fourier
mode; compare zero-forcing paths; and exercise late cancellation and retry.
These are numerical regression definitions, not experimental validation.
This entry does not add time-dependent feedback, moving-wall forcing,
prescribed reservoir forcing, stabilized complete-duration forcing, or
entropy-admitted forced composition.

## Acceleration with stationary mixture reservoirs

The native mixture and reacting owners expose
`step_with_body_force_and_reservoirs`. Characteristic reservoir admission,
caloric inversion and species/enthalpy face fluxes retain one shared owner.
Only interior cells receive mechanical half-kicks, with their actual stage
density; the prescribed exterior bath remains stationary and undepleted.
Source and reservoir checks precede trial publication, and every shortened
transport interval restarts the complete original reacting composition.

`ReactingEuler1d::step_with_dissipation_body_force_and_reservoirs` supplies
the same baths to both molecular halves and both advective stages. Mechanical
impulse/work, chemical formation-energy change, advective enthalpy and
molecular species/enthalpy/stress exchange remain distinct. Each enters the
complete momentum and absolute-energy account once. Existing final chemistry
temperature and force-Mach gates apply after the second molecular half.
No extra source workspace or new boundary/transport kernel is introduced.

`reactive/tests/forcing/reservoir.rs` covers the inert native source/face
limit, reservoir species transport, mass-derived mechanical impulse with
reaction and diffusion, zero-source parity, stationary-bath admission and
late cancellation/retry. This stationary prescribed-bath extension does not
infer finite reservoir depletion, moving exterior momentum or entropy
admission for the forced composition.

## Acceleration through complete stabilized reacting intervals

`step_with_stabilized_dissipation_and_body_force` and
`advance_stabilized_with_body_force_for/to` compose the original interior
mechanical source with chemistry, RKL2 molecular transport and optional
stationary species-bearing baths. Native source stages use actual mixture
density, and final chemistry-temperature and force-Mach admission follows
the complete molecular/chemical composition. RKL2 retains its original
nonlinear applicability and temporal-refinement requirements.

The original interval owner carries one coupling-attempt budget and one
molecular RHS budget through both halves and every rejected attempt. Accepted
mechanical impulse, conjugate work and source representation error accumulate
with compensation in `StabilizedReactiveIntervalReport::source_change`.
Only the physical impulse/work enters the original-to-final momentum and
absolute-energy balance. Representation error remains a diagnostic and
cannot become heat or conceal a conservation defect.

A zero source retains the original unforced numerical path. Zero-length calls
still admit force shape/controls, baths and the original molecular operator.
A late cancellation, source refusal, final-account failure or exhausted work
restores the entire initial conserved field and clock, including after an
accepted prefix. Existing reusable rollback/source storage suffices; the
extension retains no per-step force trace or extra integration kernel.

`reactive/tests/forcing/stabilized.rs` covers independently accelerated
reaction-diffusion Fourier transport and refinement, complete source impulse
and work, original zero-source identity, no-op admission, and whole-window
rollback/retry with stationary baths and a global work refusal.


## Reacting transport in a prescribed gravitational potential

`ReactingEuler1d::step_with_isothermal_gravity` composes the original finite-rate
chemistry half steps around `MixtureEuler1d::step_with_isothermal_gravity`.
The prepared atmosphere retains its original gas, spatial reconstruction and
constant-area reflecting geometry. Actual evolving species supply density and
calorics; reaction-induced pressure changes drive the same native momentum
update. The original Riemann species mass flux supplies gravitational work and
potential transfer. No separate rho*u work source or second gravity kick is added.

`ReactiveGravityStep` retains the native reacting report and a complete-step
potential-energy account, including both chemical halves. The gas first-law
residual includes mechanical work once; the gas-plus-potential residual restores
that work and adds the locally reduced potential change. Physical gravity
impulse and the reference-traction representation correction remain separate.
`IsothermalMixtureGravity::potential_energy_change` exposes the same compensated
local-mass-change reduction used by native transport, avoiding subtraction of
large absolute potential totals. Numerical chemical mass/energy discrepancies
remain visible; they are never repaired by adding heat.

A shortened acoustic step repeats the whole composition from accepted state.
Final caloric admission, potential representability and cancellation checks
precede publication of state, clock and report. This adds only stack accounting;
existing trial fields, reactor storage and transport scratch are reused. Zero
gravity follows the original unforced numerical path. The immutable atmosphere
is a reconstruction reference, not a constraint that suppresses reaction-driven
motion or forces chemical equilibrium. This entry does not compose molecular
gravitational diffusion, external reservoirs, moving walls, extra body forcing
or the entropy-audited lane.

Focused G0/G1/G3/G4 regressions cover hydrostatic isomer relaxation with both
gravity signs, heat-release-driven motion and native C/A/C replay, physical
formation-reference invariance, complete gas-plus-potential balance, zero-gravity
parity, source mismatch and cancellation during final gravitational accounting.


## Distributed acceleration with two inertial mixture walls

`MixtureEuler1d::step_with_inertial_walls_and_body_force` evaluates the original
forced ALE map inside every proposal of the existing simultaneous wall-speed
solve. Both force halves use their actual mixture density and physical volume;
force density and acceleration retain distinct units. Pressure impulses and
conjugate moving-wall work come from that same transport proposal and alone
supply the gas reaction on each wall. Interior body impulse/work remain explicit
external inputs in `gas.source_change` and are subtracted once from the combined
gas-plus-walls balance. Neither wall receives an implicit copy of the gas source;
its weight or other applied acceleration must be supplied as its own force.

The original acoustic/contraction bounds, force Mach checks, nonlinear solve,
per-wall work gates and combined momentum/energy admission remain authoritative.
Source refusals repeat the complete mechanical attempt. Field, geometry, clock
and both velocities publish together after the final cancellation gate. The
existing source buffer is allocated once on nonzero use and reused across
pressure proposals; there is no per-proposal field clone or second transport
implementation. Zero forcing retains the original unforced numerical path.
Values are held in moving-cell order, not interpreted as a sampled Eulerian
space/time field. This entry retains the original closed, constant-area,
adiabatic and unequal-caloric-contact limitations.

Focused G1/G3/G4 regressions cover uniform free fall with explicit gas and wall
loads, gas temperature and absent-species retention, force-density impulse from
actual moving volumes, replay through the native prescribed-wall map, separate
wall work, zero-source parity and late cancellation of all accepted states.


## Reacting gas acceleration and two finite moving walls

`ReactingEuler1d::step_with_inertial_walls_and_body_force` places the original
chemistry half steps around the forced native two-wall pressure solve. The first
chemistry half uses accepted volumes, each mechanical proposal uses the original
F/2-ALE-F/2 owner, and the second chemistry half uses actual moved volumes. Both
wall velocities are solved from finite mass and pressure reaction; gas momentum
remains resolved. A force field applies to gas in moving-cell order. Wall loads,
including each wall's weight, are separate authored inputs.

The complete gas report counts physical chemical storage, boundary work and
interior body work once. `gas_energy_change_j` restores both pressure and body
work to the source-subtracted physical energy residual. The combined gas-plus-
two-wall energy residual subtracts body work and wall-load work while retaining
exterior-pressure work; its momentum counterpart likewise subtracts body impulse
once. The native wall report still owns the middle mechanical stage only.
Source-work representation defects and chemistry discrepancies stay diagnostic,
not heat corrections or hidden wall forces.

Shorter mechanical durations restart chemistry, forcing, pressure and both wall
velocities from accepted state. The final actual source Mach gate follows the
last chemistry half, then complete momentum/energy admission and cancellation
precede atomic publication of species, mesh, clock and walls. Existing trial,
reactor and source buffers are reused. Zero forcing retains original reacting
wall behavior. This entry does not compose molecular diffusion, finite-wall heat
or a hydrostatic reconstruction reference, and inherits the original declared
reaction, splitting and unequal-caloric-contact limitations.

Focused regressions cover independent physical chemical-gas/wall energy and
momentum balances under heat release plus body forcing, analytic isomerizing
free fall with explicit wall weights, zero-source parity, and cancellation after
the candidate chemistry/geometry/wall solution has been formed.

## Reacting steady adiabatic pressure changes

The public `equilibrium_flow` owner closes a steady representative gas
stream across a passive throttle, compressor or expander using the original
element/charge-constrained ideal-gas chemistry. This extends the fixed-species
phase machinery; it does not replace its distinct real-fluid phase model.
The caller supplies original species/feed inventory, inlet temperature and
pressure, final pressure, source-valid temperature brackets, Gibbs/root work
bounds and absolute energy/entropy allowances.

The inlet is the equilibrium of the declared feed elements at the supplied
T,p. An unequilibrated feed cannot silently retain its original enthalpy through
that preparation. For machinery, the original prepared isentrope supplies the
reversible source-integrated Delta H at the final pressure. Identical Gibbs
policies and the identical source-bearing inlet are required in both SP and HP
paths. Compressor work into the fluid is Delta H_s / eta; expander signed work
is eta Delta H_s, with finite 0 < eta <= 1. A passive throttle permits only
nonincreasing pressure and declares exactly zero shaft work. A nonzero machine
pressure change must resolve work with the correct strict sign. Exactly equal
pressures return the identical inlet without an SP or HP query.

Actual products come from the existing prepared pressure/enthalpy owner.
Original NASA-9 sensible and formation energy, mixing entropy, elements,
charge, molecular-count changes and reacting derivatives remain with that
owner. Flow work is already in H; neither an extra pV term nor a second chemical
heat source is introduced. Source-relative H and S increments, rather than
subtraction of rounded absolute totals, govern acceptance. Actual Delta H minus
declared shaft work must pass the independent extensive energy allowance.
Actual chemical Delta S must be no less than minus the explicit entropy
allowance and is reported without clipping. Product mass remains checked by
the native owner against the exact original feed mass.

The accepted balance distinguishes declared shaft work, actual stream
enthalpy change and the numerical energy residual. Scaling by prescribed
positive mass flow divided by conserved feed mass yields their respective
power rates and chemical entropy production rate. Overflow and loss of a
nonzero rate to underflow are refused. This is a steady, adiabatic, negligible
kinetic/potential-energy-change model; it does not infer mass flow, geometry,
choking, finite-rate reaction, condensed phases or transport from equilibrium.

With feature `material`, `equilibrium_flow::material` consumes the existing
`EquilibriumGasEnthalpyPath` directly. Its original bank, source identity,
transport admission, conserved parcel mass and independent density-volume,
material H and composition-freezing residuals remain mandatory. The common
flow energy balance uses the actual material source-integrated Delta H.
The entropy account and reacting response are native chemical-equilibrium
quantities; frozen material heat capacity, sound speed and transport retain
their original meaning and are not substituted for reacting derivatives.

At most two original reference TP preparations and one SP plus one HP query
are made; throttling needs only the HP reference and query. Every trial retains
the original bounded Gibbs/backtracking/temperature controls. The balance
reports successful original composition and temperature iteration totals.
Cancellation is polled through the original source paths and immediately
before complete publication; borrowed models and feed remain unchanged.

The genuine-source regressions compare a single-species limit with analytic
machinery temperatures/work and stream rates, and reacting A2 <=> 2A throttles
and machinery with independent mass-action, enthalpy and entropy roots. They
also cover large formation/entropy references with small work, exact
equal-pressure identity, independent physical residual admission, original
source/budget refusals, incorrect pressure/efficiency requests and cancellation
at final publication. They are numerical conformance tests, not empirical
validation, finite-rate predictions or interval certificates.

## Distributed acceleration with molecular and thermal moving walls

`ReactingEuler1d::step_with_inertial_walls_dissipation_and_body_force` and
`step_with_inertial_walls_dissipation_heat_and_body_force` compose the existing
body-force moving-volume owner with the original simultaneous two-speed wall
root. Every proposal restarts D/2-C/2-F/2-ALE-F/2-C/2-D/2 from accepted gas;
the thermal entry wraps that sequence with the original finite-store or
prescribed-reservoir H/2 stages. Old and moved volumes remain authoritative
for their respective chemistry, molecular, heat and source stages.

Uniform or cellwise acceleration is converted with each actual source-stage
density; force density remains independent of density. The native source
receipt retains actual impulse, midpoint-velocity work and representation
residual. Only pressure and longitudinal stress supply mechanical wall
impulse/work. Direct wall weight belongs in the explicit wall loads, so gas
gravity cannot be counted twice as a wall reaction. The complete energy
balance includes physical chemical gas energy, actual finite thermal storage,
both wall kinetic changes, exterior pressure and applied wall-force work,
minus prescribed heat and distributed body work exactly once. The momentum
balance likewise removes the external body impulse once.

Shortened acoustic intervals restart all stages and the whole wall root.
Accepted gas, clock, mesh, molecular binding, both velocities and thermal stores
publish together after final chemistry, Mach, mechanical and cancellation
admission. Existing flat field/source workspaces are reused; no second
molecular, chemistry, viscous-heating or wall solver is introduced. Valid zero
source follows the original unforced arithmetic.

Focused G1/G3/G4 regressions reconstruct total energy and momentum from actual
endpoints, check native pressure/stress wall reactions, replay the original
forced moving-volume stages with a cellwise force-density/volume impulse
oracle, and exercise zero-source parity and late cancellation with finite heat
stores. These additions retain the native heat/splitting accuracy, explicit
parabolic/acoustic budgets and source-domain restrictions. They do not claim
general second order, hydrostatic well balancing, implicit FSI stability,
entropy stability or a molecular gravitational-diffusion law.

## Accelerated molecular gas in single and shared piston chambers

The original reacting molecular piston owner now accepts a held body force
through `step_with_piston_dissipation_and_body_force` and
`step_with_piston_dissipation_heat_and_body_force`. The thermal entry retains
the original finite end stores or prescribed heat reservoirs. Actual body
impulse is distinct from the native fixed support and moving-wall reaction;
physical chemical gas energy, finite thermal storage and piston kinetic
energy balance against the original boundary loads, heat and body work.

`step_with_shared_piston_dissipation_and_body_force` and
`step_with_shared_suspended_piston_dissipation_and_body_force` allow separate
uniform/cellwise acceleration or force density and Mach controls in each
original chamber. Their material banks, reference energies, molecular laws
and chemistry remain independent. The original common-speed root runs both
complete forced molecular trials at one duration, with original old/new
volumes and a single shared piston position. The mechanical balance removes
both external body impulses/works exactly once. Each fixed outer support
retains its native pressure/stress impulse. The suspension retains its
original finite spring storage and separate mechanical damper loss; no such
loss is added again as gas viscous heating.

Zero source follows the original numerical path. A shortened trial restarts
both fields and their common root; final admission and cancellation precede
publication of either field, either molecular geometry binding or the piston.
The single heated owner also commits thermal state atomically. Regressions
cover physical energy/momentum with independent force-volume and
spring/damper oracles, zero-source replay, adiabatic thermal/molecular
equivalence, and cancellation of the two-field publication. These entries
retain the original source, acoustic, parabolic, geometry and held-stage
accuracy limitations; they add no hydrostatic equilibrium, entropy theorem
or cross-piston heat inventory.



## Reacting converging-nozzle mass flow and sonic choking

`equilibrium_flow::nozzle::equilibrium_nozzle_flow` resolves pressure-driven
discharge from an equilibrium stagnation reservoir through an explicit effective
area. The original prepared chemical SP path supplies each trial's species,
density, direct source-relative enthalpy/entropy changes and reacting sound
speed. Velocity satisfies `v^2 = -2 Delta h`; the existing ordinary-branch
pressure search selects either the requested subsonic back-pressure state or
the smooth sonic limit `c_eq^2 - v^2 = 0`. Actual flow is `area * rho * v`.
Feed amounts set an elemental/charge inventory and its original reference,
not a prescribed mass-flow rate. Formation energy and flow work remain in h
exactly once. Frozen Cp/Cv or gamma never selects a reacting throat.

The private common pressure search is shared with the existing fixed-species
phase-equilibrium nozzle. Its descending scan, two-sided sonic/pressure
acceptance and sampled monotonicity checks retain the existing arithmetic.
The phase consumer retains its original PT/EOS work accounting, phase-count
admission, prepared entropy path and native equilibrium acoustics. Chemical
flow retains the original Gibbs conservation/response acceptance and bounded
temperature/backtracking solves. Neither owner substitutes another EOS.

The caller declares a connected, ordinary single-sonic branch, a positive
source-supported pressure floor, explicit physical sonic/pressure/energy
tolerances and a whole-call pressure-work ceiling of at most 256. The source
floor is a search limit; it cannot stand in for an unsupported exit or throat.
The selected sonic bracket must meet both its pressure-width and sonic-residual
allowances. Sampled monotonicity is numerical evidence, not a global uniqueness
or choking certificate for arbitrary authored caloric laws. A shallow supported
expansion does not first require the source to support an unused sonic state.
Vacuum back pressure still requires a supported positive-pressure throat.

Equal pressures or zero area return the exact original equilibrium inlet and
zero flow after one preparation. Negative area, reverse pressure drive,
unresolved positive enthalpy release, nonfinite/underflowed mass flow,
unsupported source states, insufficient pressure brackets, exhausted budgets
and cancellation refuse without publishing a partial stream or consuming a
physical reservoir. The native SP workspace is reused throughout; only inlet
and final outlet publish complete source receipts. If the final bracket selects
an earlier workspace state, an explicitly charged original query restores and
rechecks it before publication. Reported pressure/composition/temperature work
includes that query. Cancellation is checked immediately before publication.

The accepted result retains source-bound inlet/outlet states, the final chemical
SP sample, original feed mass, effective area, mass flow, speed, reacting Mach
number, choking status, source-relative stagnation-enthalpy residual and the
sonic pressure bracket. This is the steady adiabatic inviscid chemical-equilibrium
limit with negligible reservoir velocity and potential-energy change. It does
not infer a discharge coefficient, heat/friction, finite-rate reaction, shocks,
divergent geometry, recipient mixing or finite reservoir depletion. The existing
reacting rocket/finite-blowdown owner keeps its separate geometry and scope.

G1 regressions compare zero/subsonic/choked flow with the constant-Cp analytic
limit and reacting A2 <=> 2A with independent mass-action entropy, enthalpy and
acoustic derivatives. The reacting sonic station coincides with that oracle's
mass-flux maximum. G3 covers feed normalization, area scaling, choked back-pressure
independence and weak pressure drive under large energy/entropy reference datums.
G0/G4 covers exact rest, shallow source support, work/floor/input refusals and
final cancellation with deterministic retry. These are numerical conformance
tests, not empirical validation or interval certificates.



### Reacting converging-nozzle material endpoints

With `material`, `equilibrium_flow::nozzle::material` publishes the native
pressure-driven SP branch through the original
`fs-material::EquilibriumGasEnthalpyPath`. The source bank, typed inlet context,
feed elements and identical original Gibbs controls must reproduce exactly the
native equilibrium inlet and conserved feed mass. The accepted SP pressure and
source-relative enthalpy change select at most one final HP query. No alternate
chemistry, EOS, species canonicalization or transport closure is introduced.

The returned velocity comes from the **actual material** source-integrated
enthalpy decrease per conserved feed mass; actual parcel density supplies mass
flux and mass flow. The published HP chemical entropy must meet the original
SP entropy allowance. Its **reacting** sound speed supplies Mach and the sonic
residual, which must pass the original nozzle tolerance on a choked branch;
an unchoked endpoint must remain subsonic. Independent numerical bounds compare
material and native SP density and mass flux. A loose HP tolerance never
replaces any of these acceptance gates. Frozen parcel Cp, entropy audits and
sound speed remain distinct from chemical-equilibrium derivatives and entropy.

The material owner retains its original source/transport domain, species,
mass, density-volume, source-relative enthalpy and freezing admission. Equal
pressures or zero area return the exact material inlet with zero flow and no
final HP query. Native pressure evaluations retain their original budget;
material work is explicitly one reference preparation plus at most one bounded
HP query, with total actual composition and temperature iterations reported.
Final cancellation refuses publication and all borrowed sources, feed and
context remain unchanged. This is a steady ideal converging discharge, not a
reservoir-depletion step, finite-rate reaction law or inferred discharge loss.

Conformance is in `equilibrium_flow/tests/nozzle_material.rs`: reacting
subsonic/choked products, original cards/transport/context/mass, independently
rechecked source energy/chemical entropy/reacting sonic closure, original
transport-domain refusal, mismatched Gibbs controls, loose-HP rejection,
exact zero flow and final cancellation with deterministic retry. These material
feature tests require the full original dependency graph; their execution is
reported separately from the native SP/HP source-flow tests.

### Prepared reacting nozzle flow with direct local energy

The existing converging-nozzle owner now supports an internal prepared-path
consumer while keeping its public steady inlet law. Its scalar trials reuse
the original chemical isentrope and allocate source receipts only at accepted
publication. `EquilibriumNozzleFlow::specific_enthalpy_change_j_kg` is the actual
outlet-minus-current-inlet drop; `sample` may retain an earlier tank-birth
reference. The native opaque caloric anchor resolves this local drop directly.
Birth entropy, actual source residuals, reacting acoustics, pressure work
ceilings and final publication checks remain active.

Focused source tests in `equilibrium_flow/tests/nozzle_prepared.rs` cover a
weak subsonic drop after depletion that exposes rounded birth-difference loss,
reacting source/energy closure, source/normalization/stale-state refusal and
final cancellation. The shared midpoint tests independently check an exact
linear result, three-rung second-order convergence, inventory limits, exact
rest, cancellation and bounded actual-source failures. Execution of full
material reservoir integration is reported separately from these small native
source tests.

### Finite reacting converging-nozzle discharge

With `material`, `equilibrium_flow::nozzle::blowdown::equilibrium_nozzle_blowdown`
connects the existing `EquilibriumGasReservoir` to an authored fixed back
pressure and effective area over a positive physical interval. It follows
choked and subsonic discharge on the original ordinary single-sonic branch.
Every trial uses the same prepared birth isentrope, original equilibrium
response and direct current caloric anchor. It neither freezes gamma nor
re-equilibrates a newly invented reservoir reference at each rate query.

The original shared midpoint owner solves
`dm = dt*mdot(m_old-dm/2)`. Its safeguarded derivative uses
`dh_chamber/dm = c_chamber^2/m`, with `d(mdot)/dh_chamber = mdot/v_exit^2`
at fixed subsonic back pressure and `mdot/c_throat^2` at an ordinary sonic
maximum. Only the actual represented mass-rate residual accepts a root.
Initial, midpoint and endpoint nozzle regimes all undergo original admission;
only the accepted midpoint allocates a full nozzle receipt. A trial beyond
outward pressure drive supplies a zero-flux bracket bound; no physical endpoint
below ambient is accepted. Pressure crossing or excess outgoing fraction
refuses without clipping the interval, inventory or energy.

The physical endpoint and integrated exported stagnation enthalpy come from
the original finite withdrawal owner. Its source-relative `U=m*h-p*V` change
is the sole energy export; `dm*h_mid`, flow work and jet kinetic energy are
not additional debits. Elements per mass retain the original well-mixed feed
ratios while molecular composition re-equilibrates. Tank species-inventory
differences are not presented as exhaust species fluxes. Equal pressure or
zero area preserves the exact accepted tank and reports zero exchange, even
after prior depletion. Final cancellation refuses publication, and all work
inherits the original bounded composition, entropy, density, pressure-search
and mass-iteration policies with actual pressure-query counts reported.

Smooth time accuracy is second order; the midpoint mass residual and outgoing
fraction bound are not temporal error estimates. This is a rigid adiabatic
well-mixed tank discharging to an external pressure boundary, not a finite
receiver, reverse-flow connection, wall-heating law, reaction kinetics,
resolved jet or inferred discharge coefficient.

Conformance in `equilibrium_flow/tests/nozzle_blowdown.rs` covers independent
perfect-gas transient cooling/energy and three-rung time refinement, a reacting
choked-to-subsonic trajectory with elemental and total-energy closure, ambient
crossing and inventory refusal, exact rest after depletion, and final
cancellation with deterministic retry, and safe path replacement that attempts to rebase
birth entropy or loosen controls. Native execution of this original
reservoir integration is reported independently of the scalar/source tests. The public prepared reservoir rejects a replaced path with a
different birth reference, feed normalization or thermochemical policy before
any subsequent density sample or finite endpoint publication. The original
reservoir implementation lives in one dependency-light leaf and is reexported
through all existing public paths; rocket and converging consumers retain the
same physical owner.

## Prepared compact-plenum wall heat candidates

`multicomponent::network::Plenum::heat_candidate_into` applies the original
implicit gas/wall source to a stationary finite mixture inventory and an
already allocated output. It retains every canonical species density and
physical volume. The existing `WallThermal` and `WallControl` own thermal
admission, nonlinear caloric work and residual budgets; the adapter supplies
`WallDrag::None` and adds no mechanical work, reaction heat or port enthalpy.
Prescribed-temperature reservoirs and original constant-capacity finite walls
are supported. A caloric wall without its immutable material-law binding
refuses through that same source; no replacement wall EOS is inferred.

The returned `WallUpdate` carries the candidate wall temperature and original
heat, finite-wall energy and external-bath balance. A composing field must
publish wall and gas together. The last cancellation poll precedes all output
writes, and successful warm calls reuse the output species buffer. The method
changes no clock. The existing backward-Euler heat source remains first order,
including when a caller places it symmetrically around conservative transport.
Typed `multicomponent::Error::Thermal` retains original refusal details;
thermal cancellation maps directly to the existing `Cancelled` transaction.

Three focused actual-source native groups cover independent two-body implicit
heating/cooling and stiff limits, variable Cp and extreme formation-reference
invariance, prescribed-bath accounting, retained zero species and output
storage, last-poll rollback and the original extreme-conductance precision
refusal. The bounded build retains the original caloric and scalar wall
algorithms; unrelated material/shared-wall/chemistry/network orchestration and
their error payload dependencies are outside that slice. This is not a full
workspace Cargo/DSR gate or an experimental heat-transfer validation claim.

### Conservative finite receiving volumes with reacting nozzle flow

The material-feature `equilibrium_flow::nozzle::pair::EquilibriumGasPair` owns
two original `EquilibriumGasReservoir` values and one immutable chemistry/feed
binding. It advances an effective-area restriction in either pressure direction.
Both volumes may have different temperature, pressure, entropy and equilibrium
molecular composition, but must share the same original elemental ratios and
feed normalization. Mass follows the initial source EOS. This is an ideal rigid,
well-mixed, adiabatic, infinitely fast chemical-equilibrium model; it does not
mix different elemental feeds, resolve port momentum, conduct heat at equal
pressure, infer valve losses, or supply reaction rates.

The original shared midpoint mass owner solves
`dm = dt * mdot(m_d-old - dm/2, m_r-old + dm/2)`. Every trial reconstructs the
donor on its original depletion isentrope and solves the receiving volume through
the original UV owner. The nozzle therefore sees actual evolving receiver back
pressure and chemical-equilibrium sound speed. The original source-bearing
midpoint nozzle is published once. Safeguarded Newton proposals include the
receiving-pressure derivative from the same native Cv, compressibility,
expansion coefficient and sound speed; the actual mass-rate residual alone
accepts the result, with bounded bisection available.

Finite exported energy is the original donor `U_before-U_after`, including flow
work. It is never replaced by `dm*h_midpoint`. If the original reservoir reports
`C = exported_H - dm*h_d,before`, the receiver solves
`m_r,new * (u_r,new-u_r,old) = dm*(h_d-h_r+p_r*V_r/m_r) + C`, with explicit
represented-mass roundoff and its source-reference contribution retained.
Direct source-bound enthalpy differences cancel arbitrary formation datums before
arithmetic. The receiving UV target is anchored to its accepted original SP
state; the returned physical SP reservoir is independently rechecked for EOS,
energy and entropy. Receiving gas establishes a new accepted birth entropy;
a prior outflow isentrope is not imposed on it. Jet kinetic energy thermalizes
within the single carried total-enthalpy account.

Both endpoints and their one report publish atomically. An endpoint pressure
crossing, nonrepresentable positive transfer, source gap, inventory/work limit,
failed element/energy/entropy balance or cancellation refuses without clipping.
Zero area or equal pressure returns both accepted states unchanged. The physical
outgoing-fraction bound is at most one half and is not a time-error estimate.
Smooth midpoint time convergence is second order; physical time refinement
remains the caller's responsibility.

Each accepted interval also recomputes total mass, elemental inventory, energy
and entropy against the immutable initial pair using original source anchors.
The energy allowance does not grow with the number of updates. Small negative
entropy production is bounded only by the explicit numerical allowance and is
reported without clipping. The reports distinguish time quadrature from physical
energy, mass and entropy residuals.

Focused source regressions in `equilibrium_flow::tests::nozzle_pair` cover
analytic frozen receiving/cooling and swapped direction, choked time refinement,
reacting mass action and element/energy closure, formation-reference sensitivity,
exact rest, pressure crossing, work exhaustion, final cancellation, and fresh
receiving entropy plus prepared-path replacement refusal.

Six pair regressions and four native caloric-anchor regressions were authored and source-reviewed. They had not been compiled or executed at publication; focused native validation follows executor recovery.


### Source-owned primitives for mechanical trace admission

`FrozenGas::primitive_from_density_pressure` exposes the original frozen EOS
primitive at a caller-supplied density and pressure, with zero velocity. It
shares the native Riemann trace's density, temperature-domain and positive
heat-capacity admission, including its source-owned sound speed. A spatial
consumer that already holds a reconstructed velocity can therefore admit its
physical Mach number without constructing conservative energy and invoking a
caloric inversion. This adds no equation of state, approximate heat capacity,
normalization, or retained numerical state. The moving forced spatial tests in
`fs-flux` exercise the accessor through actual reconstructed-face admission.

## Prescribed-volume compact pressure-work candidates

`dynamics::network::Plenum::exchange_with_volume_candidate` extends the original
rigid extensive exchange with an explicitly supplied endpoint volume. Its first
candidate uses the original source pressure; its final SSP RK2 candidate uses
the initial and complete forward-stage source pressures with the same two-stage
port quadrature. The original thermally perfect calorics recover both pressures
and admit the resulting inventory. No constant-gamma closure, independent heat
source or additional port enthalpy is introduced.

Mass and internal energy advance extensively and are divided by the represented
endpoint volume. The mechanical receipt is `PlenumBoundaryWork`; its work uses
the represented volume difference and is positive under compression. Remote
wall work is separate from the flow work already included in transported total
energy. Equal volumes retain the original rigid numerical path. An inconsistent
forward-stage volume, exhausted inventory, unsupported caloric endpoint,
unrepresentable work or cancellation refuses without changing the original
plenum. The adapter allocates no storage and owns no clock. Its composing owner
must impose affine volume motion over the interval, source-compatible transfers,
acoustic/fractional-volume CFL and atomic publication with the connected field.

Six focused actual-source native exchange tests passed, including the retained
rigid/blowdown cases and three new groups covering an independent variable-Cv
entropy integral with second-order refinement, extreme frozen formation-gauge
invariance, separate transported and moving-wall work, and final-poll rollback.
This algebraic source does not solve piston inertia, infer wall geometry or
claim experimental validation or a complete workspace gate.


## Prepared chemical pressure feedback in compact mixture volumes

`multicomponent::network::Plenum::react_candidate_into` advances a stationary
finite inventory through the existing prepared `reactive::CellChemistry`.
The original mechanism, fixed-U,V solve, elemental/charge admission, caloric
reconstruction and formation-energy convention remain authoritative. Chemical
conversion changes the stored sensible-energy gauge by the opposite species
reference change; the adapter adds no separate heat-release source.

The method preserves volume and zero momentum, writes into an already prepared
species buffer, and returns the original `ChemistryReport` for the caller to
publish with its candidate. All output writes follow the final cancellation
poll. Successful warm calls reuse the existing reactor and output storage.
Diagnostic cell indices are supplied by the composing owner without affecting
physics. Shape, original source-domain, precision, work and cancellation
refusals leave the output unchanged.

Three actual-source native groups passed: independent reversible decay under
backward Euler and SDIRK2, exothermic temperature/pressure and chemical-energy
closure under positive/negative formation shifts, retained absent species and
prepared storage, last-poll rollback/deterministic retry, malformed output and
the original extreme-reference precision refusal. The bounded build executes
the original prepared reactor, nonlinear stages, cell conversion/caloric
reconstruction, and original fs-la decay/reaction dependency. Unrelated
network orchestration, molecular/piston wrappers and finite thermal-partner
wrappers are outside that source slice. No full workspace Cargo/DSR verdict,
coupled temporal-order guarantee or combustion validation is claimed.
