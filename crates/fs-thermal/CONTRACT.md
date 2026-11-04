# fs-thermal contract

## Crystallization in spatial heat exchange and finite molds

`enthalpy::crystallization::CrystallizationSimulation` can attach before the
first scheduled birth or adopt active all-amorphous cells from the original
`DepositionSimulation`. Bound materials receive a private
kinetic exposure; unbound molds/substrates keep their original caloric laws,
storage basis and thermal slope bounds. All cells exchange heat through the
existing link, boundary and CFL operator. Mass, geometry, cell/material IDs and
the original extensive thermal-storage account remain with that owner.
Attachment requires the original thermal law and kinetic parent branch to
carry the same complete caloric identity. Agreement at one starting state is
insufficient; unidentified or crossing-but-different caloric curves refuse.

Birth, contact, separation and boundary events use that original timeline's
admission and coverage owner, in authored order. Every thermal/kinetic interval
ends at the next event or target; target-time births activate without aging.
Each newly deposited bound cell starts with explicitly zero exposure and its
unchanged incoming enthalpy. Existing cells retain their exposure, and ordinary
molds receive no crystal history. Precrystallized inflow requires an explicit
history-transfer model and is not inferred. Empty intervals consume no thermal
steps; events do not renew the whole-advance substep budget. Incoming mass and
birth enthalpy enter their original accounts once, while heat residuals retain
cell-relative differences. A late birth-domain refusal or cancellation rolls
back the entire candidate, including its schedule cursor. Manual endpoint
boundary changes require due events to have been processed first.

Each candidate resolves crystallinity and temperature together through the
shared constitutive core. Crystal mass is part of the existing specimen mass,
and latent release is part of the existing enthalpy. Compensated audits count
both signed global residuals and cumulative absolute local caloric/roundoff
defects. Late material failures, cancellation and exhausted budgets return no
accepted history. Cloning checkpoints time, enthalpy and exposure; identical
time partitions replay bitwise. Endpoint boundary/conductance changes use the
original physical ownership and coverage checks.

`with_equilibrium_remelting` binds a prepared `EquivalentAgeRemelting` to
one original material. The full source-derived capability must match that
material's kinetic law; selection survives cloning and scheduled births.
Other materials retain their irreversible growth behavior. Existing heat
links and boundary channels can now remelt earlier crystals when later hot
material is deposited. The signed crystal-mass observation becomes negative
during net melting; it is neither an added mass source nor another heat debit.
Three timeline regressions cover latent heating followed by superheat and
recrystallization, hot incoming material remelting an older deposit through
finite contact, and original-history preservation on binding/cancellation
failure. These timeline regressions require the complete native dependency
graph to execute.

Focused regressions cover adiabatic and spatial calorimetry,
recalescence while cooling, different-Cp finite-mold exchange, timestep
refinement, source/identity/budget refusal, initial-chart roundoff, ordinary-cell
domain failures and atomic cancellation/restart, as well as independently aged
scheduled births, finite-mold contact/separation and continuation across events.
This fixed-reference-geometry,
uniform-cell reduction has first-order heat transport and explicit kinetic
timestep limits. Conductivity evolution requires a declared supplied law;
front morphology, moving geometry, mechanical unlocking and work feedback are
outside this transport model.

## Finite-rate crystallization and latent thermal feedback

`crystallization::CrystallizationLaw` retains an independent Nakamura exposure
`z = integral K(T) dt` and absolute crystal mass fraction
`X = Xmax * (1 - exp(-z^n))`. Its two original caloric branches share one datum:
`h = (1-X) h_parent(T) + X h_crystal(T)`. Distinct heat capacities and the
original temperature-dependent enthalpy gap enter this one account; latent heat
is never added again as an independent source. Relative conversion and maximum
absolute mass crystallinity are separate parameters.

Isothermal stepping returns the required thermostat heat. Nonisothermal
stepping integrates exposure with explicit midpoint while reconstructing both
midpoint and endpoint temperatures on the same enthalpy account. Work and
caloric tolerances are bounded, cancellation publishes no candidate, and
out-of-domain or unsupported remelting trajectories refuse. Immutable source owners may
supply positive whole-domain heat-capacity lower bounds to accelerate caloric
inversion; generic laws retain safeguarded bracket convergence. Sixteen focused
native core Rust tests pass, covering analytic JMAK, adiabatic recalescence with
unequal heat capacities, second-order time refinement, accelerated affine and
nonlinear inversion, datum shifts, endpoint behavior and cancellation.

The declared source domain may extend above melting for completely amorphous
feed. Both original caloric branches and the kinetic source must cover that
full domain; no chart is extrapolated. The supplied Nakamura rate must be
exactly zero at and above melting, and positive queried rates refuse instead
of being clamped. Once any nonzero crystalline mass fraction exists, direct
state queries and enthalpy inversion restrict temperature to the original
domain at or below melting. Cooling superheated amorphous feed enters the same
conservative crystallization account. The default growth path refuses
remelting; the separate capability below explicitly selects that model. The core regressions include superheated
cooling, source-domain/rate refusal and crystal-preserving melting-boundary
checks. A growing-timeline regression deposits 550 K feed against a finite mold
and follows sensible cooling into crystallization without a second latent heat
source; it requires the full native thermal dependency graph to execute.

The default is an irreversible isokinetic constitutive model with explicit
caller parameters and applicability. It supplies no material calibration, stochastic
nucleation, grain/dendrite morphology, segregation, remelting/reset or mechanical
locking. Time integration error remains distinct from caloric closure.

### Explicit equilibrium remelting and thermal cycling

`CrystallizationLaw::equivalent_age_remelting` admits and caches the original
melting-temperature calorics in an immutable `EquivalentAgeRemelting`.
Both branches and the original zero-at-melting kinetic source must cover that
temperature. `advance_specific_enthalpy_with_remelting` preserves the original
enthalpy at its midpoint and endpoint while latent heat is absorbed on that
plateau. Surviving crystal MASS fraction is obtained from the original forward
Avrami relation after a bounded equivalent-age reconstruction; complete melting
returns zero effective age and the original amorphous branch. Cooling resumes
the original growth law from that age.

This separately chosen history approximation supplies no retained-nucleus,
lamellar-population, dissolution-rate, molecular melt-memory or mechanical
unlocking model. Smooth unmelted paths retain the original midpoint method;
transition steps require refinement and carry no second-order event-accuracy
claim. A projected midpoint that would suppress continuing undercooled growth
refuses. Zero heat and zero growth rate preserve an accepted plateau state
bitwise, including an admitted caloric residual, without another inverse solve.
Five native core regressions cover partial/complete melt, superheat and renewed
growth, conserved energy, source/domain/history refusal, large datums and
saturated ages, false-arrest prevention, cancellation and stationary holds.

Layer L3. Shared conservative thermal transport and solid/vapor exchange.
Runtime dependencies: std, fs-la, fs-math, fs-solver. Spatial discretizations and material
constitutive data remain owned by their existing crates.

## Drifting Maxwellian incident molecular moments

`sublimation::molecular::drift::incident_maxwellian` supplies incoming mass,
conditional vector velocity and crossing energy relative to the original
species enthalpy at a locally moving planar surface. Normal and tangential
drift use the same truncated three-dimensional Maxwellian. The caller supplies
partial density, original specific gas constant, gas temperature, relative
bulk velocity and an already-unit outward normal. Empty species yield exactly
zero mass and no conditional moments; unrepresentable nonzero outputs refuse.

Central moments use deterministic Gaussian functions and compensated arithmetic.
Away-flow tails use positive integral-ratio recurrences with checked terminal
bounds and at most 496 backward iterations. Scaled logarithmic amplitudes retain
representable low-probability incidence; dimensional full moments handle cold
inward beams once the omitted tail is below floating-point resolution. There
is no physical Mach cutoff. The signed mass change from zero normal drift is
evaluated independently, preserving small net deposition near saturation even
when the gross beams round equally. Energy includes tangential slip and the
flux-weighted normal second moment; adding bulk kinetic energy to a stationary
beam would not give these moments.

The kernel adds no state, thermochemical datum, film resistance, continuum
stress, velocity-dependent sticking or Knudsen-layer model. Original material
owners must combine these incident moments with emission/reflection and the
laboratory momentum/energy transformation. Focused tests specify independent
quadrature, opposite-half reconstruction, rotations, tiny drift, scaled tails,
cold beams, vacuum and representability refusals.

## Conducting ambient faces and separate physical baths

`ConductiveAmbient` and `ConductionObserver::with_ambient_faces` bind a selected
cell's complete original exposed ambient area to the same finite half-cell
material and declared surface-resistance law as conducting reservoirs:
`G(h) = A / (d/k(h) + R_surface)`. The original ambient and reservoir baths
remain separate and may change temperature through boundary events. Geometry,
reference mass (except explicitly bound planar phase layers), material identity
and coverage stay fixed; changing a bound
face refuses. Zero conductivity insulates exactly. Zero surface resistance
prescribes temperature at the physical face while retaining the material
half-cell resistance.

The explicit observer protocol carries total ambient conductance through the
original boundary heat and CFL owner. Authored htc is preserved, ambient heat
is not reclassified as an applied source, and the entropy owner retains the
exact applied coefficient and actual bath temperature. Observer composition
rejects duplicate channel ownership. Existing advection, flow, electrical,
rendering and solid-feedback wrappers forward this channel to their original
inner observer.

Implicit chains and graphs use the same boundary-aware residual and analytic
conductivity derivative for each bath, including nonzero dk/dh on a latent
plateau with zero dT/dh. Existing unselected channels retain their original
behavior. Source-only implicit observers explicitly refuse ambient overrides
rather than silently dropping physical heat. Source, topology, work-budget
and cancellation failures publish no changed thermal history.

## State-dependent contact in implicit thermal networks

`enthalpy::conduction::implicit::contact::ImplicitContactLaw` supplies positive
area-specific contact resistance and its two specific-enthalpy derivatives on
the original ordered caloric objects. Chain and graph `with_contacts` attachment
requires unique active links and zero authored resistance on those links, so
source and fixed contact resistance cannot be counted twice. Unselected links
retain the original fixed-resistance path. Explicit zero resistance requires
zero tangents; isolation still checks the source domain.

The existing endpoint supplies `G=A/(d_a/k_a+R''+d_b/k_b)` to the same conservative
backward-Euler residual. Each Jacobian column includes both its bulk-k response
and `-G²/A * dR''/dh_i`, then divides by that cell's original mass for the native
energy unknown. The resistance-fraction form avoids explicitly squaring G.
Opposite heat rows retain opposite signs. No additional heat source, enthalpy
storage, clock or nonlinear solver is introduced. The existing nonsymmetric
graph solve, backtracking, entropy, energy gates and transactional refusal own
acceptance. Interface geometry, pressure and interfacial storage remain fixed.

## Existing thermal timeline

The enthalpy, advection, contact, observation and event APIs retain their
documented control-volume semantics. A cell stores total energy as mass times
specific enthalpy. Every material chart in one balance uses a compatible energy
datum. The shared timeline is re-exported by fs-conduction; consumers must not
implement independent copies of its storage update.

`advance_observed_continuing_limited` continues the accepted observer from an
earlier `advance_observed_limited` interval within one enclosing physical
advance. It uses the same thermal kernel and all original observer hooks,
without calling `begin` again. Coupling owners can split heat and mass exchange
without resetting a source's per-window energy or physical history. The caller
must supply the remaining shared substep allowance and retain the original
deep-clone/rollback contract; continuation does not renew that work budget.

`DepositionSimulation::apply_heat_j` adds signed endpoint heat to those same
cells without advancing time or re-integrating sources and contacts. The
original material law resolves `h_new = h_old + Q/m`, including latent storage.
Same-time scheduled events must be processed first. The method checks actual
cell-relative stored energy against the supplied joules, with separate local
and aggregate absolute-error tolerances; opposite lost increments cannot hide
behind zero net heat. Compensated sums preserve small changes beside large
caloric reference energies. Domain, arithmetic and cancellation failures return
no candidate. The energy audit identifies the input as external to the thermal
subsystem; a coupling owner must account for the corresponding upstream debit.
No mass, geometry, boundary, event cursor or accepted time changes, and no
mechanical-to-heat conversion fraction is inferred. Five focused regressions
cover latent/signed inputs, unchanged clocks/sources, large reference datums,
local/aggregate loss gates and transactional refusal.

### Enthalpy transported by balanced mass flow

The existing upwind graph evaluates transport relative to one canonical
start-state enthalpy per connected flow component. Fixed mass makes the common
reference cancel analytically; its numerical value cannot create heating or
erase a small source by subtraction of gross formation-energy powers. Component
identities are prepared once, and the same relative fluxes feed thermal and
optional entropy balances. Latent enthalpy still moves at equal temperature.
The gross absolute transported-enthalpy diagnostic retains its original datum.

The existing cell mass-balance admission remains 64 machine epsilons of the
incoming/outgoing rate scale. Inside that allowance the equations enforce zero
continuity; the time-integrated absolute cell mass residual and its
`|h_reference * mass_residual|` energy equivalent are reported explicitly as
roundoff diagnostics. Neither is a physical heat source, mass loss, or second
energy debit. G1/G3 regressions cover a flowing calorimeter at large positive
and negative reference offsets, small heating in a uniform circulating loop,
disconnected latent-heat loops, canonical edge orientation, admitted mass
roundoff, and matching energy/entropy reference invariance. This does not add a
pressure, velocity, variable-mass or chemical-reaction closure.

`Advected::with_selective_outlets` explicitly binds a carried composition and
specific enthalpy to existing outlet IDs before stationary species admission.
The original boundary owns the mass rate. On acyclic graphs the same species
solver subtracts those component fluxes, checks availability without clamping,
and propagates the remaining mixed stream. Every selected outlet declares the
same caloric reference as the registered cell and inlet laws. Its carried phase
enthalpy is not inverted through the remaining gas law. The existing thermal
source, continuity, CFL and external-enthalpy accounts include that outlet once;
no separate latent-heat source is invented. Replacing boundary or selective
data invalidates the previous stationary species solution.

The constitutive caller owns selection rates, phase equilibrium and source
alignment. No condensation kinetics or new thermodynamic data are inferred by
this transport primitive. Negative remaining species flux, absent remaining
mixed flow, an unknown/duplicate boundary, a mismatched caloric reference,
selective recirculation and the optional entropy audit refuse. Focused G0/G1
regressions cover component conservation, downstream propagation, phase
enthalpy outside the gas inversion domain, the original thermal energy account,
and these unsupported conditions. This does not add transient species storage.

### Liquid enthalpy through a stationary solid skeleton

`Advected::with_liquid_carrier` retains the original total-cell `mass*h`
storage while a locally balanced mass graph carries only the liquid phase's
enthalpy. For an admitted pure-substance isothermal latent plateau, liquid
carries the liquid-side plateau enthalpy during coexistence and its original
specific enthalpy when fully liquid. The component reference is evaluated on
that same carried-enthalpy function. Circulation at the melting temperature
therefore preserves stationary solid fractions, while incoming superheat can
melt material through the existing energy account without a second latent
heat source or an independent phase-fraction integrator.

The descriptor requires the complete original phase chart and retained
caloric identity, affine liquid mass fraction on the plateau, equal constant
phase densities and fixed cell mass consistent with density times volume.
The public generic constructor requires the caller to provide the original
complete chart; selected temperature samples do not establish that binding.
Every flowing donor and receiver must contain liquid. Local mass
balance, conservative pair transfers, the existing combined explicit CFL and
implicit energy residual remain authoritative. The carried-enthalpy slope is
zero on the plateau and one in the liquid branch, so the original total-mass
outflow bound remains conservative.

Original `MassBoundary` inlets and outlets can supply balanced open liquid
streams. Every inlet uses that same material and enthalpy datum, lies within
the original source chart and is entirely liquid at or above its liquidus
enthalpy. Entry and discharge cells require liquid support and the same
fixed-density inventory. Outlets carry the current liquid enthalpy, including
in partially solid cells. The original external-enthalpy balance counts the
net stream heat once. Isothermal flushing preserves the stationary solid;
inlet superheat fills its latent-energy deficit before sensible heating.
Boundary and carrier attachment order does not change admission or arithmetic.

The native `MaterialFlowFeedback::with_liquid_carrier` adapter passes the
original pointer-bound phase curve's full knots and retains the descriptor
through explicit graph rebuilds and implicit Stokes endpoints. It requires
the existing physical creeping-flow heating owner. Native velocity represents
liquid volume flux per total face area; homogeneous retained inertia refuses
in both attachment orders. Open ordinary mobile-mushy flow uses the existing
opposed pressure-mouth work owner and explicit pure-liquid enthalpy on actual
inflow, including backflow. Its fixed thermodynamic reference pressure remains
separate from the hydrodynamic pressure gauge: adding hydrodynamic `p/rho` to
the inlet enthalpy would count the physical pressure supply twice. Exact-solid
pressure-mouth topology remains unsupported. Selective streams, stationary species,
bulk-enthalpy entropy auditing and nonempty bulk-enthalpy reconstruction also
refuse until their own phase-specific closures are supplied. Alloy melting
ranges, unequal phase densities and moving solid skeletons are outside this
model.

Five thermal regressions cover isothermal circulation, large datum shifts,
implicit latent melting by superheat, fully liquid equivalence and admission
failures. Two native flow regressions exercise repeated physical wall-heating
epochs and the implicit Stokes energy equations. Five additional open-stream
thermal regressions cover phase-specific discharge, reversed flow, large
enthalpy datums, latent-to-sensible timestep refinement and refusal. Two
pressure-driven flow regressions retain analytic Darcy heat and melt-dependent
resistance across explicit rebuilds and implicit endpoints, including pressure
gauge shifts. These integration tests
require the full native thermal/flow dependency graph to execute.

## Sublimation and deposition

### Deposited coatings on curved substrates

`phase_geometry::RetreatGeometry::{spherical_coating,cylindrical_coating}`
extend the existing mass-owned front with a fixed inert substrate. The phase
mass and `FrontGeometry::volume_m3` describe only the coating; `extent_m` is
its thickness and the exchange area is the outer shell area. Cylinder end
caps remain inert. `occupied_volume_m3` additionally includes the substrate
and is the volume needed for dilute-body/gas-displacement admission.

The geometry follows the shell identities
`V = 4*pi*t*(r*r + r*rc + rc*rc)/3` and
`V = pi*L*t*(2*rc+t)`, where `r=rc+t`. Inversion and recession use factored
forms so a thin deposit is not lost by subtracting nearly equal radii or
total volumes. Each query uses constant storage and bounded scalar work;
there is no measured hardware-speed claim.

The original thermostatted front integrates coating thickness at
`dt_thickness/dt = -j/rho`, refreshes area from the same conserved mass, and
reports extinction exactly at the substrate. After extinction, phase area
is zero and unused time is returned to the caller. This does not remove the
substrate or assign it a replacement boundary law. Existing finite-energy
moving-chamber consumers can use the same geometry without a second front
integrator.

G0/G1/G3 tests in `phase_geometry` cover analytical shell mass/area/tangents,
thicknesses below one ULP of the substrate radius, deposition heat and mass
closure, subdivision invariance and coating-only extinction. These are
uniform constant-density, shape-preserving coatings with an existing seed;
they do not infer nucleation, porosity, substrate thermal storage, curvature
corrections to equilibrium pressure, or gas-volume work.

`sublimation::SublimationLaw` supplies a bounded immutable saturation-pressure,
latent-enthalpy, density and molecular-mass chart. The material crate supplies
reference ice/CO2 correlations and an explicitly parameterized
Clausius–Clapeyron approximation. A numerical correlation is not a sourced
specimen certificate, and supplied density/accommodation/film coefficients
retain their caller-declared applicability.

Positive mass flux means solid becomes vapor; negative flux deposits vapor.
The input pressure is the volatile partial pressure. Outgoing molecules use
surface temperature and incoming molecules use gas temperature. One common
accommodation coefficient preserves detailed balance at equal temperature.
The optional dilute gas film is combined with interface kinetics by series
resistance. Zero accommodation or film coefficient seals the surface.

`surface_flux` returns signed mass and latent-energy fluxes, their analytic
surface-temperature derivatives, and outward-normal interface velocity.
`isothermal_mass_transfer` integrates an ideal, thermostated finite vapor
reservoir exponentially at fixed surface temperature. It uses the shared
deterministic expm1 implementation to preserve short timesteps and does not
overshoot its frozen-temperature equilibrium. Its caller supplies the thermal
work needed to hold the surface at that temperature.

`advance_surface` couples that transfer to the caller's original T(h) chart
with a bounded implicit endpoint solve. The endpoint vapor carries h_end + L.
The energy equation is:

    m_initial * (h_end - h_initial) + transferred_mass * L(T_end) = 0

This is algebraically the balance of remaining solid energy and transported
vapor enthalpy. The latent account is part of vapor energy, not another
independent energy debit. All mass limits are applied before deriving heat,
vapor inventory and recession. A shifted common caloric datum leaves the
physical temperature and mass solution unchanged.

The nonlinear solver requires a continuous nondecreasing caloric chart and a
bracketable endpoint in the declared constitutive domain. It refuses
unsupported temperatures and failed convergence rather than clamping them.
The default budget permits 96 evaluations; callers may select 1 through 4096.
Each evaluation polls the supplied checkpoint. Functions return candidate
values and mutate no caller state on failure.

### Deposition on an initially bare thermal substrate

`sublimation::coating::advance_supported_surface` separately conserves the
substrate and deposited-solid inventories while sharing one interface
temperature. The deposit may start at exactly zero mass when the caller
explicitly admits heterogeneous deposition. Only the deposit can sublime;
exhaustion leaves the substrate intact and permits later regrowth. The
substrate's original T(h) chart and deposit's original h(T) chart remain
authoritative, with independent enthalpy reference levels. The bounded solve
closes `M_s*(h_s1-h_s0) + m_d0*(h_d1-h_d0) + dm*L(T1) = 0`, where positive
dm enters the vapor. The vapor carries the deposit's `h_d1+L`; substrate,
deposit and vapor energy changes close one common account.

This is a thin, isothermal deposit approximation with frozen exchange area
and the existing prescribed-temperature reservoir semantics. It supplies no
nucleation barrier, porosity, resolved coating temperature gradient or
contact-resistance model. The coating timeline below applies this transaction
with the original heat-flow owner. Other spatial consumers must explicitly own
their geometry and contacts; whole-cell sublimation remains a separate model.
Six focused tests cover seed-free growth, finite-vapor inventory and energy,
extinction/regrowth, inactive boundaries, independent caloric reference shifts,
invalid charts/inventories, cancellation and bounded nonlinear work.

### Supported coatings on the native thermal timeline

`sublimation::coating::heat_supported_surface` applies signed heat at fixed
substrate/deposit masses by solving
`M_s*(h_s1-h_s0) + m_d*(h_d(T1)-h_d(T0)) = Q`.
It reuses the bounded endpoint solver, including backtracking when the
substrate-only outer trial leaves the domain but the combined-storage solution
is admissible. Independent material datums remain intact. Unrepresentable heat
is either explicitly bounded by the declared residual tolerance or refused.

`enthalpy::sublimation::coating::CoatingSimulation` owns this storage together
with the original thermal timeline and supported phase operator. Substrate
cells retain their original material, mass, enthalpy, geometry and identity;
deposited mass and its enthalpy remain separate, including exact zero inventory.
The wrapper's mass/energy queries include BOTH condensed inventories. Each
selected surface explicitly supplies its original pure-solid h(T), compatible
phase law, kinetics and independent prescribed/thermostatted vapor reservoir.
Future births may carry an authored coating; its mass and sensible energy enter
the incoming account exactly once.

The material adapter `fs_material::thermal::PureSolidThermalLaw::
specific_enthalpy_at_temperature` supplies original solid h(T) directly from
its borrowed curve. It retains exact source knots, varying sensible slopes and
the melting-onset endpoint, with no constant-Cp surrogate or extrapolation.
Ambiguous all-solid isothermal states and unrepresentable inverse queries
refuse; mixed/liquid states never extend the solid branch's support.

Native heat links, ambient and reservoir exchange, source power, contacts,
separation and boundary events retain their original owner. Heat storage
includes the coating before the following phase update on every CFL substep.
The substrate's existing CFL bound remains conservative for continuous
nondecreasing deposit h(T). Endpoint coverage changes cannot retroactively
alter deposition, and future cells have no pre-birth exposure. Deposit
exhaustion preserves the substrate and its subsequent ability to collect vapor.

The complete candidate includes event time, substrate enthalpy, deposit mass
and vapor inventory. Refusal or cancellation publishes none of them. One
full-advance thermal budget and bounded local nonlinear solves apply. Local
energy differences and compensated sums preserve caloric reference handling;
actual mass-rounding corrections are charged to the energy residual, and
opposing local defects cannot cancel out of its absolute admission.

This remains first-order heat/phase splitting with fixed reference area and
conductances and a shared substrate/coating temperature. It does not resolve
thick-layer resistance, porosity, melting, moving geometry, dynamic
observer-based sources, or a shared gas chamber. Finite baths in different
surface bindings are independent. Seven native integration tests cover combined
storage/conduction/source work, zero-inventory growth and regrowth, finite vapor
conservation, exact event exposure, coated births, cancellation/budgets,
restart and independent caloric datums; seven heat tests additionally cover
nonlinear storage, domain backtracking and unrepresentable heat.

## Numerical and physical limits

The coupled endpoint is a first-order reduced model. Area remains frozen over
one phase step; recession uses endpoint solid density. Finite gas temperature
is prescribed: its thermostat owns subsequent sensible-heat equilibration.
The dilute film does not model non-dilute Stefan flow, pore resistance, plumes,
curvature, nucleation, melt competition, radiation or unresolved gradients.
Those mechanisms require explicit compatible models, not hidden coefficients.

For a whole-solid body, returned zero mass denotes extinction. Its spatial
consumer must update topology, contacts and conservative remapping or refuse
that endpoint. A supported coating can instead reach zero inventory while
retaining the substrate and its topology, as described above.
Reporting recession does not itself move a mesh, alter a level set or change a
flow/electrical geometry map. Uniform-cell use requires a justified internal
temperature approximation. Refining time cannot repair an invalid spatial
model.

The scalar flux and finite-reservoir kernels allocate no memory and perform
constant work. The implicit solver uses constant memory and at most its
declared evaluation budget. This is an algorithmic statement, not a measured
hardware throughput or cross-ISA reproducibility claim. Transcendentals use
the existing fs-math owner; current executed evidence is same-ISA.

## Focused verification

G0/G1/G3 tests cover signed deposition, equilibrium, unequal gas temperature,
analytic tangents, kinetic/film limits, finite-reservoir exact relaxation and
semigroup behavior, tiny timesteps, inventory limits, implicit cooling and
heating, mass/energy/recession consistency, caloric-datum shifts, first-order
refinement, cancellation and exhausted nonlinear budgets. Material tests
compare the actual correlations to published reference values and their
analytic derivatives to finite differences. These tests establish the
declared reduced model; they are not experimental validation of a complete
sublimating object or a general spatial Stefan solution.


## Phase exchange on the existing heat timeline

`enthalpy::sublimation` composes the surface operator with the original
`DepositionSimulation`. `sublimate_cell` atomically changes the active cell's
actual mass and specific enthalpy plus its vapor inventory, without advancing
the clock. `advance_sublimating_to` performs bounded first-order heat/phase
splitting with one total thermal work budget and one global mass/energy audit.
It honors birth and boundary/contact times exactly, uses the exposed area from
the start of each interval, and never sublimates a not-yet-born cell. The
existing caloric chart remains authoritative; rejected out-of-domain trial
points do not invalidate a successfully converged in-domain endpoint. All
failures and cancellation preserve the original timeline and reservoir.

The whole cell must remain a pure solid on its original caloric branch.
Temperature alone cannot establish solid fraction on a melting plateau;
mixed-phase inventories require a separately owned partition. The supplied
finite reservoir belongs to this one surface and remains thermostated. A
joint chamber shared by several surfaces needs a joint backpressure solve.
Reference geometry/contact areas remain fixed and recession is reported.
Complete cell exhaustion refuses until a consumer owns topology removal.
There is no implied moving mesh, gas-energy solve or material calibration.

Vapor energy includes the solid's original sensible-enthalpy datum plus latent
heat; the separately reported latent account is never an additional debit.
Ten new native integration regressions cover balances, frost growth, covered
surfaces, birth/contact timing, nonlinear chart trials, transactional refusal,
whole-advance budgets and first-order refinement. The combined exact-source
thermal/material suite passed 95 tests with one existing performance test
ignored. `cargo run -p fs-material --example sublimating_ice` demonstrates the
real ice law with explicit illustrative calorimetry in a maintained vacuum and
a finite chamber; both native example cases passed mass/energy assertions.


## Several surfaces sharing one vapor chamber

`sublimation::chamber::advance_chamber` solves one shared endpoint vapor
pressure together with the original conservative surface enthalpy operators.
Its backward-Euler mass balance admits warm-to-cold redistribution, including
condensation of vapor emitted during the same interval into an initially empty
chamber. Depositing surfaces never receive independent copies of the original
vapor inventory. One volatile species with constant molar mass, pure-solid
caloric branches and compatible enthalpy references are required. Gas volume
and temperature stay prescribed.

A frozen-temperature predictor and safeguarded scalar pressure interpolation
avoid a dense coupled matrix; each trial costs O(surfaces), with explicit
surface, local-iteration and total pressure-trial ceilings. Unsupported outer
trial bounds are backed away from. The same accepted transfers reconstruct
final vapor mass and drive every solid mass, latent/sensible-energy and
recession account. The nonlinear pressure/mass closure residual is reported
separately from actual published mass conservation. Refining time controls the
first-order integration error; the residual does not certify model accuracy.

Eight new native regressions check independent simultaneous algebra, transient
refinement, stationary-pressure circulation, initially empty chambers, shared
inventory, area partition/order, refusal and finite-capacity energy closure.
The combined thermal/material suite passed 113 tests with one existing benchmark
ignored. Thermally evolving gas, molecular heat accommodation, mixed phases and
mesh/topology changes remain separate physical owners.

### Relative entropy on the shared thermal timeline

ThermalLaw may provide its exact relative integral of dh/T at the original
fixed-pressure/composition convention. Missing entropy remains unavailable;
no substitute heat capacity or entropy reference is invented. The opt-in
advance_to_with_entropy and BackwardEulerStep::evaluate_with_entropy use the
original accepted thermal states, actual interface coefficients and the same
start/endpoint flux quadrature. EntropyAudit separates exact material storage,
actual ambient/support reservoir entropy, conductive production, unclassified
source energy divided by cell temperature, and the numerical balance defect.
The source flag survives cancellation of signed sources. Birth entropy is
excluded on both sides; subsequent storage starts from actual birth enthalpy.
A trial audit is accepted only with its trial. Ordinary energy-only APIs do no
entropy integration. G1/G3 tests cover isolated exchange, refinement, reservoirs,
latent storage, equilibrium, deferred birth, resolved implicit sources, missing
material entropy and cancellation without publication.


## Shared chambers on actual multi-cell heat histories

`enthalpy::sublimation::chamber` applies one joint phase exchange to canonical
cell identities using each cell's original caloric law and exposed area.
`sublimate_chamber` validates every candidate before committing any solid or
vapor state. `advance_chamber_sublimating_to` composes original heat transport
with joint exchange, honors start-of-interval exposure and exact scheduled
birth/contact boundaries, and retains one whole-advance thermal work budget.
The report checks combined solid/vapor mass and the original total energy
account, while keeping each surface's recession separate.

Whole cells must remain pure solid; gas temperature and reference geometry
remain fixed. Complete extinction requires a separate topology owner and
refuses atomically. Ten native integration regressions passed for inventory,
redistribution, original caloric charts and heat graph, canonical order,
event timing, work ceilings and late cancellation. The runnable fs-material
example `ice_cold_trap` uses the real ice law to compute frost growth on a cold
film supplied by a warmer film through an initially empty shared chamber.


## Molecular counterflux for nonisothermal interfaces

`sublimation::molecular` returns both nonnegative half-space mass fluxes and
one signed net flux from the original bounded constitutive and kinetic owners.
The equal-temperature path retains stable pressure-difference evaluation.
Different surface and local vapor temperatures remain distinct; volatile
partial pressure is never replaced by total carrier pressure. Scalar work is
constant and allocation-free. Four native regressions cover equilibrium,
unequal temperatures, pressure semantics and refusal boundaries.

An energy consumer must retain the opposing streams: zero net mass does not
imply zero heat exchange. For Maxwellian ideal vapor with equilibrium internal
modes, each stream carries h(T)-R_specific*T/2 in an explicitly aligned energy
reference. The material gas adapter owns that caloric evaluation and both
endpoint energy accounts. This interface assumes energy/angle-independent
equal sticking and specular nonsticking reflection; it supplies no gas-film,
diffuse thermal accommodation, Knudsen-layer or moving-geometry model.


### Molecular normal momentum and explicit Maxwell reflection

`sublimation::molecular::momentum` supplies allocation-free half-space momentum
moments for the SAME admitted gross mass streams. The crossing-molecule mean
normal speed is `sqrt(pi*R_specific*T/2)`; emission and sticking loads add even
when net phase transfer is zero. Reflected nonsticking molecules are a distinct
account, including carrier species. Explicit Maxwell diffuse fractions retain
incident mass and reemit that fraction at the surface temperature; the rest
reflects specularly. A thermal accommodation measurement alone does not select
an angular scattering model. At equal temperatures the complete load recovers
the original gas pressure independently of the sticking probability.

Finite input, probability, overflow and nonzero-underflow admission applies.
These moments add neither heat nor a second ambient-pressure load. Their
stationary surface and zero-drift local Maxwellian assumptions do not supply
moving-wall, continuum-stress or gas bulk-momentum evolution. Six regressions
cover independent Gaussian quadrature, equilibrium and sealed limits, vacuum
recoil, zero-net counterflow, Maxwell reflection and extreme-scale refusal.


## Original thermal storage for external solid/gas coupling

`enthalpy::sublimation::external` exposes the original pending events, work
budget and cell caloric evaluation, and transactionally applies prepared solid
mass/enthalpy endpoints to the existing thermal storage. Every full expected
cell, positive remaining inventory, caller-admitted pure-solid h interval,
original caloric endpoint and original-geometry exposure is checked before a
single write. Endpoint contact events may change current exposure; they cannot
retroactively alter the preceding interval's declared area. A surrounding
coupled owner mutates its private trial without another graph clone per phase
step. Joint gas/solid conservation and physical phase admission remain with
that owner; this boundary does not invent calorimetry or geometry evolution.

The material gas timeline now composes this boundary with the actual original
heat operator and one jointly evolving gas. Four native owner tests cover
state preservation, stale/bad batches, endpoint covering and late cancellation.
Both existing isothermal sublimation timelines also jump directly to the next
event/target while the WHOLE graph is empty, avoiding unbounded tiny pre-birth
intervals. Two regressions prove a long pre-birth gap with zero thermal work
allowance. The native sublimation suite passed all 51 selected tests.

## Heat-driven condensed-fuel pyrolysis

`enthalpy::pyrolysis::PyrolysisSimulation` extends the original activated
thermal cells with separate virgin and char inventories. It retains their
identity, geometry, contacts, boundary conditions and extensive storage, and
uses the same conductive/boundary heat operator through a private caloric
callback. No second conduction graph or alternate transport algebra is used.
`CondensedCaloric` adds source-owned h(T) to the original `ThermalLaw::T(h)`;
admission checks that the charts agree. Virgin and char may have distinct
nonlinear sensible calorics. `VolatileEnthalpy` supplies the escaping gas's
physical enthalpy on a compatible formation reference.

One first-order non-oxidative reaction converts virgin material into an explicit
char mass yield and outgoing volatiles. Its accepted endpoint temperature
determines the backward-Euler inventory, retained char and outgoing enthalpy:

```text
k(T) = A exp(-E/(R T))
m_virgin_new = m_virgin_old / (1 + dt k(T_new))
extent = m_virgin_old - m_virgin_new
mass_out = (1 - char_yield) extent
energy_out = mass_out h_volatile(T_new)
```

The local energy equation balances the change in remaining solid storage plus
that exact outgoing enthalpy against supplied heat. Reaction enthalpy is already
in the compatible material references; it is never a second source or debit.
The independent virgin/char caloric changes and reaction enthalpy differences
are evaluated before gross formation energies are combined. Per-cell releases
retain their integrated mass, enthalpy and consumed virgin mass for downstream
gas coupling, even when emission temperatures change within one advance.

Tiny rate-duration products use stable loss fractions. Published consumption
is reconciled to actual retained inventory, so an unrepresentable decrement
cannot emit fictitious mass or heat. The absolute kinetic mass-roundoff and
reaction-energy equivalents are reported separately and checked against the
caller budgets; they are numerical discrepancies, not physical transfer streams.
Local, cumulative absolute and whole-advance mass/energy residuals cannot be
hidden by cancellation between cells. Errors and cancellation publish no
candidate. Boundary and conductance changes use original transactional admission.

This initial lane requires endothermic single-step kinetics, monotone sensible
caloric charts, explicit temperature/source domains and fully activated fixed
topology. It refuses pending geometry events, exhausted work budgets and
unrepresentable depletion. It does not seed residual solid mass. Conduction
retains the original explicit stability bound using the current condensed mass
and mixture capacity; local reaction is implicit in its final temperature.
Near depletion this fixed conductive mesh can become stiff. Moving surfaces,
porosity/permeability, pore-gas transport, oxidation, latent-phase state,
multi-step wood chemistry and fire-spread closure remain separate capabilities.
The coupled method is first order, with no unconditional stability or entropy
theorem. Conductivity changes require an explicit caller-supplied law.

Eleven focused regressions cover independent isothermal and adiabatic
solutions, distinct char calorics, first-order temperature/mass refinement,
equivalence to original heat transport with reaction disabled, spatial heat-
driven mass release, caloric-datum shifts, small-rate representation, and
domain/cancellation/work-budget refusal. The native full thermal library suite
passes 164 tests with one existing performance test ignored. This is focused
native Cargo evidence, not a repository DSR gate or experimental validation.

`heated_pyrolysis_slab` exercises the actual 1D graph with a resolved front
half-cell heat balance, irradiation, emission and convection. It uses retained
NIST/MaCFP scalar PMMA inputs in an authored heating scenario; its optional
20% char sensitivity is explicitly not a wood material fit. `macfp_tga` uses
the same kinetic update against the retained MaCFP FDS material-temperature
history. Source data, units, calibration/verification distinctions and commands
are in `data/reference/compartment_fire/pyrolysis/README.md`.

### Implicit conduction and finite-rate decomposition

`pyrolysis::implicit` removes the explicit conduction stability limit while
retaining the original physical model. Its nonlinear unknown is absorbed cell
energy. Every trial calls the existing `solve_cell` for temperature, virgin/char
inventory and escaping volatile enthalpy, then the original
`heat_rates_with_calorics` evaluates endpoint conductive/Robin/source powers.
The actual residual is absorbed energy minus duration times endpoint power.
No alternate heat graph, heat capacity, kinetic law or heat-of-reaction source
is introduced.

`implicit_workspace` prepares original dense LU and numerical buffers for at
most 128 cells. The Jacobian combines the original link/boundary conductances
with local caloric-response probes, avoiding a complete graph reevaluation for
every column. Original bounded heat-rate maps and candidate/result objects
still allocate; this is not an allocation-free whole-advance claim.
`advance_implicit_to` returns the unchanged original release/audit format plus
actual nonlinear, local-solve, linear-solve and backtracking work counts. Local
and whole-advance conservation gates include original kinetic roundoff and the
separately reported graph heat-transfer summation discrepancy. Cancellation or
nonlinear/domain refusal publishes no candidate.

Seven focused native tests cover an independent stiff two-capacity endpoint,
simultaneous heating/conduction/reaction, char and caloric-datum invariance,
actual Robin coverage, temporal refinement, activated Arrhenius comparison
with the original fine explicit method, and refusal/cancellation/reuse. The
stiff case completes in one implicit interval while the explicit owner exhausts
its four-substep budget. This remains first-order backward Euler with caller
time-resolution controls; removing the conduction CFL is not an unconditional
nonlinear-convergence, temporal-accuracy or entropy theorem. Geometry and
authored coefficients stay fixed within each original advance.

### Conservative noncharring slab recession

`pyrolysis::slab::RecedingSlab::new` adopts the original homogeneous x-aligned cell
chain with caller-sourced constant density, conductivity and cross-section.
Actual retained mass sets thickness through `L = m / (rho A)` with a fixed
back face. This is the single-material density/length relation in the pinned
FDS Solid Chapter's shrinking/swelling section; fixed-back anchoring and
conservative rezoning are explicit choices of this implementation.

Each bounded interval uses the original implicit conduction/reaction owner.
Excessive actual per-cell mass loss retries the complete unpublished interval,
without accumulating rejected volatile releases. Accepted specific enthalpy is
reconstructed with limited linear slopes in retained-mass coordinate and
integrated over overlaps onto equal-mass cells. The original caloric inverse
recovers temperature; original cell geometry and series conductances are rebuilt.
Mass and extensive enthalpy are transferred, not averaged temperatures. Cell IDs
are spatial slots after rezoning rather than tracked material particles.

Original volatile mass/enthalpy integrals are unchanged and retain their
production slots and homogeneous source binding. The original aggregate audit
includes actual remap storage/discrepancy once; separate mass/energy residuals,
rounded-geometry mass equivalent and enthalpy redistribution are returned.
The redistribution diagnostic is not lost heat or a PDE error bound. No seed
mass, heat-capacity floor or release clipping extends the physical domain.
Geometry and coefficient evolution remain first-order split from each implicit
interval. The charring extension below supplies separate residue geometry.
Pores, mechanical work, thermal expansion, gas-domain volume displacement and
complete burnout require separate models.

Seven native groups verify independent isothermal inventory/recession,
nonlinear-caloric overlap conservation, monotone reconstruction/refinement,
formation-datum shifts, reaction-disabled identity, complete retries and
cancellation/replay/source/work refusals. The combined native thermal suite
passes 178 tests with one existing performance test ignored; DSR is unavailable
in this execution environment, so this is not a repository-level gate result.

`heated_pyrolysis_slab` retains its original five arguments and fixed default;
an optional sixth `receding` argument enables the new noncharring path with the
actual front half-cell resistance. Original NIST/MaCFP scalar inputs and the
effective volatile enthalpy remain unchanged. The 120-second 24/48-cell runs
complete, with emitted-mass spatial sensitivity of about 1.88%; their actual
accounts are retained in
`data/reference/compartment_fire/pyrolysis/receding_slab_example_results.json`.
These authored heating cases do not establish experimental PMMA validation.

### Conservative charring slab shrinkage and swelling

`RecedingSlab::with_char` adds retained residue to the same finite-rate
virgin-to-char-and-volatile material, implicit heat solver and release accounts.
`CharProperties` supplies independent char density, conductivity, source and
a declared `ConductivityMixing::Parallel` or `Series` arrangement. The
existing `new` constructor remains the noncharring path.

The declared condensed volume is `V = m_virgin/rho_virgin + m_char/rho_char`.
The back remains fixed; front motion may be positive recession or negative
recession (swelling). This additive-volume choice is explicit, not the general
conditional FDS shrinking/swelling algorithm. Parallel conductivity is the
volume-fraction arithmetic mixture used by the pinned FDS Solid Chapter;
series conductivity is its harmonic counterpart for the separately declared
unresolved arrangement. All rebuilt internal links and the public
`cell_conductivity_w_m_k` and `front_half_cell_conductance_w_k` queries use the
same current component mixture.

Equal-retained-mass rezoning conserves total mass, char mass and physical
enthalpy, with virgin mass recovered as the complementary component. Char
fraction and sensible enthalpy use a paired limited reconstruction that
preserves donor means, nonnegative component inventories and the common
caloric-domain bounds. Formation references travel with their component
masses. The original scalar energy owner recovers mixed temperature with
reaction disabled during inversion; rezoning creates no volatile stream or
additional reaction heat. Separate virgin/char remap residuals prevent total
mass cancellation from hiding composition drift.

Frozen-geometry intervals are bounded by both the existing mass-loss limit
and `maximum_cell_volume_change_fraction` (default 0.25). This also controls
swelling with little mass loss. An excessive change retries the entire
unpublished interval, retaining no rejected release. Reports include maximum
volume change, component and total mass discrepancy, physical-energy
discrepancy, geometry roundoff and original caloric evaluation work.
Cancellation or failed budgets return no candidate.

Five additional native groups establish independent component inventory,
volume, fixed-back motion, outgoing enthalpy and both conductivity limits;
nonlinear-caloric remap and series-link consistency; better-than-threefold
remap refinement with independent phase energy-datum shifts; swelling retries;
and source refusal, cancellation and deterministic reuse. All 12 slab tests
pass, and the focused native thermal library passes 183 tests with one existing
performance test ignored. This is not a full-workspace DSR result.

Geometry/coefficient evolution remains split and first order in time, and
remap refinement does not certify global PDE accuracy or an entropy theorem.
No gas-volume displacement or pressure work is inferred. Pore transport,
permeability/porosity evolution, char oxidation, cracking, multi-step wood
chemistry, latent phases, complete burnout transitions and an experimentally
validated wood material remain outside this closure.


## Geometric conduction with evolving material conductivity

`enthalpy::conduction::ConductionObserver` prepares existing thermal links with
explicit centre-to-face distances, contact resistances and immutable scalar
`ConductivityLaw` bindings. It queries each participating cell once per actual
explicit substep and shares that value across incident faces. Internal links
reuse the original series resistance owner. A surface reservoir adds the
half-cell material resistance and an explicit areal film/contact resistance;
zero external resistance prescribes temperature at the face, not at the cell
centre. Zero conductivity explicitly insulates. Geometry, reference mass by default,
material identities and face bindings remain fixed and are checked on reuse.
Material selection and mixed-phase microstructure remain the constitutive
owner's responsibility.

The original observer protocol now accepts step-local reservoir conductances.
Both pairwise links and reservoirs enter the same original heat-rate and CFL
calculation. Their authored values are restored before accepted-step hooks and
events. No reservoir heat is reclassified as an applied source. Composed
observers reject conflicting reservoir ownership, and the existing advection
wrapper forwards the hook. `advance_observed_with_entropy` forwards each
physical hook once and retains the exact applied sources, links and reservoirs
for the original entropy account, including actual bath temperatures. General
internal transfers and applied sources retain explicitly unclassified entropy
contributions; temporal entropy defects remain distinct from physical
conduction production. Failure or cancellation publishes no partial physical
state.

The prepared graph uses O(cells + faces) storage and
O((cells + faces) log(cells)) work per substep with ordered maps. This is an
algorithmic bound, not a measured throughput claim. The existing explicit
coefficient-freezing and caloric CFL applies; it is not an implicit solver or
a time/spatial error certificate. Focused G0/G1/G3 tests cover series/Robin
resistance, dynamic boundary CFL and entropy, exact event splitting, independent
nonlinear transient and Stefan front/heat refinement, composition, insulation,
domain/work limits and transactional cancellation.

### Mass-dependent resistance of planar phase layers

`ConductionObserver::with_planar_phase_layers` explicitly binds selected cells
to constant bulk density and a fixed planar footprint. Current conserved mass
sets each normal half-thickness to `d = (m / rho) / (2 A)`. The original
Fourier/contact resistance owner uses this distance for both link sides,
reservoirs and ambient faces, while the original caloric and CFL owner uses
the same current mass. Growing frost therefore changes both storage and
resistance; its conductivity still comes from one original query per cell.
Unselected cells retain the fixed-mass refusal.

Admission checks the reference mass against density times reference volume,
the total two-face area against `2 A`, authored normal half-distances and
per-side contact measures. Source identity, reference mesh, contact coverage
and footprint remain fixed. A nonzero positive mass/thickness must remain
representable. All selected channels must be normal to those planar faces;
this is an explicit caller assertion, not an inferred geometric orientation.
The fixed-mass implicit chain/graph constructor refuses this explicit mode.

This closure supplies no mesh motion, sidewall growth, porosity, nucleation,
complete extinction, gas displacement work or density evolution. In particular,
mixed-phase calorics do not justify constant mixture density automatically.
The existing first-order thermal/phase split requires timestep refinement.
The four G0/G1/G3/G4 groups in `conduction/phase_layer_tests.rs` exercise growth
and retreat resistance, both link sides, finite heat storage/conservation,
temperature-dependent conductivity shared by both baths, per-side coverage,
unchanged default admission, cancellation and the implicit-mode refusal.


## Implicit enthalpy transport on geometric conduction chains

`enthalpy::backward_euler::BackwardEulerStep` additionally accepts optional
endpoint reservoir-conductance replacements through
`evaluate_with_reservoirs` and `evaluate_with_entropy_and_reservoirs`.
Original bath temperatures, areas, heat rates and entropy remain authoritative;
authored coefficients remain unchanged through evaluation and acceptance. Existing
coefficient callbacks retain their original behavior. Unknown cells, absent
baths, negative/nonfinite conductances and cancellation refuse transactionally.

`enthalpy::conduction::implicit::ImplicitConductionChain` solves fixed-mass,
nearest-neighbor chains using that same residual, not a second heat-rate owner.
Every caloric and conductivity law is bound to the original borrowed object.
A fused constitutive query supplies T(h), k(h), their analytic local derivatives
and an explicit admitted enthalpy interval. Physical latent plateaus retain
zero dT/dh while their conductivity may continue changing with phase inventory.
The exact endpoint Jacobian differentiates both the temperature difference and
the half-cell/contact/surface conductance, in cell ENERGY coordinates.

Only three matrix diagonals are assembled. The reusable `fs-la` adjacent-pivot
workspace keeps the linear solve O(cells) in work and storage; original
map-based graph/constitutive bookkeeping remains. Newton candidates whose
nonzero increment cannot change the stored enthalpy are represented by exactly
zero increment before evaluation. This does not bypass or relax the original
underflow policy, heat rates, per-cell residual or total-energy acceptance.
Every resulting candidate is reevaluated by the original owner.

Declared controls bound timestep, accepted intervals, Newton iterations,
backtracks, interval halvings and total residual evaluations, including
rejected trials and final entropy queries. Bounded backtracking respects the
material domains; nonconvergence or a singular linear system may trigger a
whole-interval retry. Scheduled events split intervals at their exact times.
Only fully accepted intervals contribute heat or entropy. A complete call
also checks independently summed original stored energy against all supplied
heat. Failure may alter scratch but publishes no physical history. The report
retains original heat/entropy, maximum accepted local residual and work counts.
Every admitted caloric law must provide the original entropy integral.

This path removes the explicit diffusion CFL restriction for an accepted
implicit interval. Backward Euler remains first order; interval halving here
is nonlinear globalization, not a time-error estimator. Global uniqueness,
unconditional nonlinear convergence, a general graph solve, evolving geometry
and constitutive accuracy beyond the supplied laws are not claimed.

Five focused chain regressions cover stiff exact solutions, latent traversal,
original-residual Jacobian differences, nonlinear time refinement, event
splitting, cancellation/budget rollback, singular falling-conductivity retries,
and a 96-cell diffusion tail with sub-ULP energy updates. Two added reservoir
regressions check original heat and actual-bath entropy plus refusal/rollback.
The full pinned-main native thermal source passes 191 tests, with one existing
ignored test; required workspace/DSR gates remain a separate unavailable check.

## Additional material conductances on the native explicit timeline

`DepositionSimulation::limit_timestep_for_additional_conductances` composes
nonnegative extra W/K with an already admitted explicit step using the original
cell masses and registered global caloric dT/de bounds. It adds the largest
extra decay rate harmonically to the current rate allowance, rather than taking
two independent CFL minima. An empty or zero map leaves the proposed step
bit-for-bit unchanged; zero caloric slope adds no restriction. Unknown cells,
invalid rates, unrepresentable steps and checkpoint failures refuse.

This helper applies no heat and changes no thermal state or authored boundary.
Call it before event clipping through the existing observer timestep hook.
It can be conservative when different cells or a caller maximum step dominate.
It bounds frozen heat-transfer coefficients, not their nonlinear derivatives
or time-discretization error. Focused native regressions cover combined
conductive/cooling limits, sequential composition, zero slopes and refusal.
The airflow material observer uses the exact hA from its start-state source
evaluation, so changing fluid data changes both cooling and native subcycling.

## Shared endpoint ambient and reservoir conductances

`backward_euler::{BoundaryConductances, BoundaryStepCoefficients}` pass endpoint
total ambient/reservoir G directly into the original heat-rate and entropy
owners. The boundary-aware residual and external-power APIs retain prescribed
temperatures, reference coverage, authored h/G and cell watts. Replacement G
is applied once; no boundary input is mutated or reclassified as generated heat.
The same G participates in signed physical watts, outgoing-conductance/CFL
accounting and boundary entropy Q/T. Existing ordinary callbacks use empty maps;
the reservoir-only APIs delegate their maps through this same general ingress.
Finite thermal stores remain original internal links, not prescribed baths.

Unknown cells, absent reservoirs, negative/nonfinite G and nonzero Robin G
without exposed reference coverage refuse. Three added native regressions cover
independent finite-cell energy/entropy/CFL, unscaled applied watts and empty-map
parity, invalid coverage and final cancellation/replay. These tests are compile
checked in this change; their Cargo runtime build was killed while compiling
the existing asupersync dependency. The earlier native test counts above are
prior validation, not runtime evidence for this extension. Geometry-model
validity and nonlinear/time accuracy remain the composing caller's responsibility.

## Declared thermal storage and closed-volume admission

`ThermalLaw::storage_basis` declares whether the original caloric law stores
specific enthalpy or fixed-volume specific internal energy. Existing laws retain
the enthalpy default. The inherited `specific_h_j_kg` field carries the selected
coordinate; the timeline applies the same conservative heat exchange and
`integral de/T` entropy account without converting one thermodynamic potential
into another. Each law can admit its own mass/geometry relationship through
`validate_cell`, called before a birth changes inventory.

`material_storage_basis` exposes the original registered law's declaration.
Nonzero prescribed mass flow and inlet calorics require enthalpy storage. The
single-surface, joint-chamber and external-endpoint sublimation paths also require
enthalpy before constructing or publishing open-mass exchange. Batch admission
finishes before any endpoint mutation. Zero-flow connections remain available
for heat-only internal-energy storage.

Focused regressions exercise two rigid bodies exchanging heat with conserved
internal energy and increasing thermal entropy; invalid mass/volume births;
forward/reverse circulation and boundary-flow refusal; and all three sublimation
entry paths. Mixed enthalpy/internal-energy batches preserve the original solids,
vapor, mass, stored energy and clock on refusal. This adds no open-mass internal-
energy model, pressure-work conversion or equation of state. Such conversions
remain the explicit responsibility of the physical material/equipment owner.

## Adaptive implicit thermal accuracy

`enthalpy::conduction::implicit::adaptive` adds first-order step-doubling
control around the original `ImplicitConductionChain`. One full backward-Euler
interval is compared with two actual half intervals; only the conservative
half-step trajectory is accepted. Comparison trials never contribute physical
heat, entropy or accepted time, and no Richardson-extrapolated state is formed.
Unequal representable half spans receive their first-order error correction.

`enthalpy::conduction::implicit::adaptive::AdaptiveConductionStepper` supplies
a typed endpoint seam for the existing adaptive controller. The original
conduction chain and material-integrated layered endpoints share the same
step-doubling comparisons, change-based storage/temperature/phase error scales,
exact event grid, original energy/entropy accounts and whole-call work ceilings.
An adapter must retain original cell/material identities, charge every residual
evaluation and completed nonlinear iteration including refused trials, and
return accepted original backward-Euler endpoints without internal physical
subdivision. This interface adds no conductivity, caloric law, extrapolated
state, global-error certificate or trajectory adjoint.

`ImplicitConductionLaw::time_observables` supplies original temperature and an
optional original phase mass fraction. The per-cell estimator checks stored
specific energy, temperature and, when requested, phase inventory. Relative
scales use changes from the accepted state, rather than absolute enthalpy or
Kelvin values. Phase accuracy refuses when an original observable is unavailable
rather than inferring phase from temperature. The h-labelled coordinate retains
the declared thermal storage basis.

Controls retain the original local/whole-call joule gates and bound comparison
spans, comparison count, accepted half steps, consecutive reductions, Newton
work and the shared residual-evaluation budget. Scheduled events are hit exactly.
Source refusal, exhausted work, cancellation and unresolved accuracy return no
partial trajectory. Reports separate local time-error estimates from algebraic
energy residuals, count all attempted nonlinear work, and expose the next proposed
comparison span for continuation.

Four native regressions cover analytic cooling and caloric-datum invariance;
phase-driven refinement on an exactly isothermal latent plateau; non-dyadic bath
events and fine-only heat/entropy; and atomic refusal for missing observables,
minimum-step limits, trial/residual budgets and late cancellation. A local
step-doubling indicator is not a certified global time error or spatial error
bound. The existing fixed-chain geometry and material-domain restrictions remain.

## Implicit conduction on branched thermal graphs

`enthalpy::conduction::implicit::ImplicitConductionGraph` extends the same
fixed-mass enthalpy solve to branched, cyclic and disconnected graphs, including
parallel interfaces and asymmetric covered contact areas. It delegates to the
original backward-Euler residual, nonlinear line search, reservoir/event owner,
local and global energy gates, entropy audit and adaptive step-doubling logic.
Cells, caloric laws, conductivity laws and geometric bindings remain original.

The graph Jacobian retains the analytic derivatives of both T(h) and k(h),
including zero dT/dh on latent plateaus and energy-coordinate mass factors.
It stores one diagonal and two off-diagonal contributions per physical link.
The existing `fs-solver` FGMRES owner applies this generally nonsymmetric
operator with diagonal preconditioning; no dense global matrix is assembled.
Each operator action is O(cells + links), and bounded restart storage is
O(cells * restart + links). Original map-based physical bookkeeping remains.
The direct `ImplicitConductionChain` path retains its O(cells) tridiagonal solve.

`ImplicitConductionLinearControls` bounds restart dimension, restart cycles
per nonlinear direction and total Krylov iterations for the whole advance.
Coarse, fine, rejected nonlinear and rejected adaptive trials share that
iteration ceiling; `linear_iterations` reports completed work. Linear
nonconvergence may use the same interval-halving path as nonlinear refusal.
Cancellation and exhausted budgets publish no physical candidate. Linear
residual convergence never replaces the original joule acceptance gates.

Regression coverage includes stiff star equilibration with unequal capacities,
independent endpoint elimination, energy and entropy accounts, latent branches,
active cycles and parallel links, an isolated cell, asymmetric covered areas,
the original-residual Jacobian and transpose, adaptive analytical diffusion,
shared work budgets and cancellation. Geometry and reference masses remain
fixed; backward Euler remains first order. No unconditional nonlinear
convergence, spatial accuracy, hardware throughput or experimental validation
beyond the supplied material laws is claimed.

## Evaporation on the original thermal timeline

`enthalpy::evaporation` connects the continuum liquid/gas-film owner to the
actual `DepositionSimulation` cell mass and original caloric enthalpy.
`evaporate_cell` changes these two coordinates atomically and leaves the clock,
identities, reference geometry, contacts and authored sources intact.
`advance_evaporating_to` composes original heat transport followed by implicit
phase transfer on event-clipped intervals. The interval's initial exposed area
survives endpoint covering, and no phase exchange occurs before cell birth.
The whole request shares the original thermal work limit.

`LiquidCellExchange` and `apply_liquid_exchanges` admit externally solved
finite-atmosphere endpoints on those same cells. Full expected-cell equality,
unique identities, original caloric/storage and cell constraints, admitted
pure-liquid enthalpy range, positive remaining inventory and interval area are
checked for the complete batch before its final cancellation gate and writes.
The external physical owner remains responsible for joint liquid/gas balances.

The caller admits the solver's specific-enthalpy bracket as a pure-liquid
branch. Original births and every thermal substep are checked against that
branch; an internal-energy law refuses open-mass evaporation. The original
caloric law supplies all temperatures, including during the phase solve.
Signed transfer permits condensation without clamping humidity. Actual mass
change, full carried vapor enthalpy and thermal storage close one whole-window
account. Latent energy is separately visible only as a subaccount of vapor
enthalpy; it is never added as a second thermal source. Numerical balances use
the unchanged absolute whole-window joule allowance.

G0 unit regressions cover latent cooling and condensation warming, actual mass
and energy conservation with conduction and heat sources, event exposure,
birth timing, original branch/storage refusal, cumulative work limits and late
cancellation. A G1 regression refines the complete heat/phase split, which is
first order. These checks do not establish experimental model accuracy.
The atmosphere is prescribed: its temperature, pressure and vapor fraction do
not evolve here. Finite atmosphere energy/species evolution requires its own
owner accepting the complete transfer. Geometry is fixed reference geometry;
recession is reported, with no inferred VOF, liquid-flow, dryout, boiling,
freezing or contact/topology transition.

## Caloric identity for constitutive laws on native stored energy

`ThermalLaw::caloric_identity` optionally exposes the immutable source-law
identity, including its specific-energy datum, as 32 bytes without adding a
hashing or material dependency to this thermal kernel. `None` means that the
law has not supplied this binding. `DepositionSimulation::material_caloric_identity`
queries the actual registered owner and refuses unknown material IDs.

Constitutive consumers indexed by native specific energy can require this
identity together with the existing material and storage-basis checks. The
identity supplies no new caloric model and does not authenticate empirical
data. Existing thermal laws and temperature-only consumers remain valid when
the optional identity is absent; consumers requiring it must refuse absence.

## Physical source transitions on the observed native timeline

`DepositionObserver::next_transition_time_s` optionally supplies an absolute
physical source/state transition. The explicit native owner queries it after
same-time authored events, before preparing an active-cell interval. Returned
times must be finite and strictly future, even when beyond the requested target.
The native interval ends at the earliest event, transition or target after all
conduction and observer timestep restrictions. A later harmonic decay limiter
therefore cannot turn a source edge into an asymptotic sequence of steps.

Transitions are not authored events and do not consume their budget. They retain
the existing substep budget and transactional candidate semantics. Coupled
observers validate each child's result before selecting the minimum; generic
thermal, optical, flow and mechanical wrappers forward the hook. Default None
preserves existing stepping arithmetic, and the existing empty-grid skip is
unchanged. Focused observer regressions check exact coupled boundaries, invalid
time rejection with rollback, coincident native events and split continuation.

## Observed absorption sources in implicit conduction

`ImplicitConductionChain::advance_observed_to` and the corresponding graph
entry point accept original `DepositionObserver` source owners on the same
native timeline. The existing backward-Euler residual, analytic conduction/
contact/cooling Jacobian, Newton line search, chain direct solver and graph
Krylov solver remain authoritative. There is no second thermal integrator.
Original source-free nonadaptive advancement shares this event/retry path.

Start-state source watts and final-interval source watts enter the original
residual exactly once, alongside authored boundary heat. The source sees the
accepted interval START state; its supplied watts remain fixed during Newton.
Conduction and native boundary cooling are evaluated at the endpoint. This is
an IMEX discretization for temperature/enthalpy-dependent absorption, not a
fully endpoint-coupled optical solve. Optical temporal refinement remains
necessary even when implicit diffusion permits a much larger thermal step.

Native events and source transitions cut intervals exactly. Observer timestep
restrictions remain active. A nonlinear interval retry re-evaluates the source
for its actual shorter duration, so a linear ramp retains the correct mean.
Source `begin` runs once per advance; work used by rejected preparations stays
charged. Preparation may replace pending sources/caches but must not publish
physical history. The source's acceptance hook runs only after the original
thermal endpoint passes its cell and total joule gates. The entire thermal and
optical candidate remains transactional if a later source check refuses.

This first entry point supports external sources only. Nonempty observer link
conductances, reservoir conductances or conservative transfers refuse explicitly;
their physical effects are never silently dropped or duplicated. Authored
boundaries are restored before native solving and acceptance. Native external
heat already contains observed absorption, so component optical accounts are
parallel descriptions of the same energy, not additional thermal credits.
Generated-source entropy retains its original unclassified status. The existing
adaptive step-doubling path remains source-free.

G1/G3/G4 checks cover independent stiff chain/graph forced solutions, implicit
Robin cooling with authored/start/interval sources, nonlinear halving of ramp
power, cumulative rejected-work limits, exact event/source cuts, split restart,
source timestep bounds, unsupported-hook refusal and late cancellation/history
rollback. Cross-crate tests additionally use actual conductor surface pulses,
state-dependent n/k response and Beer-Lambert volume deposition. Those checks
retain independent native and optical energy accounts; they do not certify
spatial, constitutive or nonlinear-optical accuracy from a small energy residual.

## Limited axial enthalpy reconstruction on native transport

`enthalpy::advection::reconstruction::LimitedEnthalpyReconstruction` prepares
bounded three-cell stencils on contiguous original centreline geometry. Minmod
limits physical one-sided enthalpy derivatives using actual adjacent lengths;
the signed native mass flow selects the donor's appropriate face. Each face
trace lies between its donor centre and neighboring sample, and an extremum
has zero slope. Missing stencils preserve the original donor-cell trace.

`Advected<()>::endpoint_transport_with_reconstruction` first runs the original
endpoint admission and source calculation, then corrects each selected internal
transfer with the same power on both cells. Original component-datum handling,
fixed mass continuity, caloric laws, external stream enthalpy and backward-Euler
energy accounting remain authoritative. Prepared mass, geometry and material
bindings are checked on every candidate. Topological junction gradients,
external streams on stencil centres, and entropy/species/selective-stream
closures refuse instead of silently borrowing an unsupported reconstruction.

The plan bounds authored stencils and retains only native cell identities and
prepared geometric factors; no new physical timeline or accepted-state clamp is
introduced. The caller solves the nonlinear endpoint residual and retains the
original cell and total joule acceptance gates. Regression sources exercise
nonuniform affine traces, smooth periodic face-error refinement, reversal and
datum invariance, an independent implicit cell balance above explicit CFL, sharp
extrema, invalid ownership/topology and cancellation. These checks concern
spatial face reconstruction and stated discrete equations, not a blanket TVD,
global convergence or arbitrary-source maximum-principle claim.

