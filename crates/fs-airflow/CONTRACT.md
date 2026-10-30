# CONTRACT: fs-airflow

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
Existing solid IQN, parallel-pump and recirculation entry points retain their
original explicit omission budgets; they do not silently acquire this model.

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
an existing forced-convection card or the Churchill-Chu vertical-plate card to
the original finite-cell enthalpy owner. Both laws share one source, exposed-area
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

`medium::convection::deposition::source::evaluate_source_cell_convection`
re-evaluates the original fluid owner at an explicitly selected reservoir or
mean film temperature for every proposed thermal endpoint. The property
temperature never replaces the actual reservoir temperature in buoyancy or
signed heat flow. Density, heat capacity, conductivity, viscosity and known
signed isobaric expansion share one checked source state. Natural convection
requires sourced expansion; forced convection does not invent it when absent.

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

`source::observer::SourceConvectionObserver` supplies these heat rates to the
original explicit `DepositionSimulation::advance_observed` timeline and composes
with other `DepositionObserver` owners. Every actual substep reads its starting
wall/film state and exposed area; accepted endpoint samples retain freshly
resolved material owners. Scheduled events invalidate the sample, and the next
substep uses the updated contact coverage and boundaries. The observer's signed
interval heat is a share of the existing external-heat audit, never a second
energy addition. Endpoint source refusal or cancellation discards both proposed
histories. Resolver cloning must preserve the native observer's rollback rules.

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
common-pressure admission bound; only an axial series path identifies that sum
with its end-to-end drop. No pump or flow-feasibility claim is added. The new
network regressions independently reconstruct cell and external enthalpy
balances, check recirculating species, exercise material/pressure substitutions,
and compare the original axial API with the same graph and wall laws.

`medium::convection::finite_gas::ForcedGasSurface::isochoric_body_candidate`
connects the same material/correlation evaluator to the generic finite-gas,
uniform-body thermal callback in `fs-couple`. Each coupled trial resolves h from
the actual candidate gas endpoint; accepted transport and correlation sources
remain beside the same gas state used for energy and acoustic observations.
The original body solver retains conductivity, radiation, Biot admission and
enthalpy, and the single exchanged heat cancels in the combined energy balance.
A positive patch must cover the whole body surface; zero area explicitly
insulates it. Bulk-gas property evaluation and prescribed relative flow retain
their correlation limits. Dittus-Boelter selects its heating/cooling branch
from the proposed signed heat inside the same residual. Its neutral zero-heat
seed uses the smaller branch coefficient, preventing an authored direction
from falsely rejecting the physical branch at the Biot gate. Constant-heat-flux
duct correlations refuse on this isothermal-body path.
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

### Source-preserving film properties in finite gas/body coupling

`finite_gas::film::FilmGasSurface` evaluates forced or natural convection at
the mean of the actual bulk gas and wall temperatures, while the existing
finite-gas/body/piston residual owns both heat debits, gas pressure work and
mechanical motion. The film state supplies properties only; the original bulk
gas retains mass, energy, composition and acoustic state. Rayleigh forcing and
signed exchange retain the full wall-to-bulk temperature difference.

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

### Natural convection in finite gas and coupled surface networks

`medium::convection::finite_gas::NaturalGasSurface` evaluates the shared
vertical-plate correlation from the coupling owner's actual gas and body query
states. The retained ideal-gas EOS supplies beta = 1/T; the same resolved
thermochemical gas supplies density, Cp, viscosity and conductivity. Pressure
and composition therefore propagate through Ra, Pr, h and the original signed
gas/body heat residual. The existing piston, phase-enthalpy and gas owners
advance their own energy, including latent heat and piston work.

`finite_gas::network::exchange::GasHeatLaw::NaturalConvection` carries that
same law into the material electrical/contact/radiation/thermal network and its
chamber-flow step. Natural, forced and explicit conductive laws may occupy
disjoint exposed patches, with aggregate area admission and one energy owner.
Accepted samples retain the exact queried gas and wall states. Candidate solves,
atomic flow-step publication and cancellation use the existing network path.

This is an unconfined isothermal vertical-plate approximation evaluated at bulk
gas properties; it does not resolve buoyant velocity or supply enclosure
corrections. The caller must admit the geometry and temperature contrast.
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
permeability supplies only the unchanged tensor/series sum of thickness over normal
permeability divided by flow area. The existing hydraulic owner evaluates current
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
feedback and parallel redistribution/mixing.
