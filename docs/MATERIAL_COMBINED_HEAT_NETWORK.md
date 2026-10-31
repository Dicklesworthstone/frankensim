# One material assembly for electrical and thermal exchange

`fs_couple::circuit_enthalpy::contacts::heat::MaterialHeatNetwork` composes an
existing `ElectricalContactNetwork`, optional `MaterialContactLink` records,
and an optional `MaterialRadiation` enclosure. It introduces no material law,
Kirchhoff solver, radiosity solver or time integrator.

```rust,ignore
let assembly = MaterialHeatNetwork {
    electrical: &electrical_contacts,
    thermal_contacts: &thermal_contacts,
    radiation: Some(&enclosure_binding),
};
let selected = assembly.prepare(&cx, &thermal)?;
let proposal = selected.step_candidate(
    &cx, &thermal, target_s, &circuit, &source_ports, controls, &energy_guess_j,
)?;
let (next_thermal, next_circuit, report) = proposal.into_parts();
```

Every evaluation uses the same trial enthalpies. Fresh bulk conductances receive
thermal contact resistance in series exactly once. Signed radiative powers are
added exactly once. The existing circuit owner deposits bulk and electrical
contact losses into their original cells. Authored independent heat sources,
convection, reservoirs, event handling and latent storage stay active in the
original thermal owner. Neither a rejected trial nor a source query advances
accepted charge, flux, cell enthalpy or the event cursor.

Preparation retains the original bulk, electrical-interface, thermal-interface
and emissivity claim selections. It does not cache physical properties. Each
subsequent evaluation queries those claims at the actual temperature and phase.
An unsupported selected source refuses rather than switching to another model.
Preparing again is an explicit new source-selection operation.

Radiating cells must use their existing bulk cards, caloric laws and complete
physical conditions. A radiation binding cannot quietly give the same body a
different pressure, temperature-axis meaning or liquid-fraction convention.
The original radiation and contact owners retain geometry, source and surface
coverage admission. All original source receipts and per-mechanism reports are
available in `report.source.material`; the complete bulk material snapshots
remain in its `material` member for downstream consumers.

## Directional heat conduction

For an orthotropic material, declare `DirectionalConductionLink` records from
`circuit_enthalpy::material` and use
`assembly.prepare_with_directional_conduction(&cx, &thermal, &directions)`.
Each original link side supplies three principal conductivity property keys,
their explicit selection policy, a nonzero material-frame identity, and a unit
temperature-gradient direction expressed in those material axes. All links
touching one cell must agree on its principal roles, selected claims and frame;
their directions can differ. Link endpoint order, area and half-cell distances
remain with the original thermal geometry.

The projected coefficient is `k_n(T) = sum_i n_i^2 k_i(T)`. Both sides enter the
existing series-resistance law, followed by any declared thermal contact.
Principal claims are selected once, retained through source continuation, and
queried at every candidate cell state alongside its electrical properties.
The normal material sample carries all three principal receipts and the applied
link coefficients. An orthotropic-only cell can set its legacy scalar
`MaterialCell::conductivity` binding to `None`; no substitute scalar is invented.
Changing orientation or principal conductivity data therefore changes heat,
temperature-dependent resistance and the original coupled electrical history.

This also composes with electrical interfaces and enclosure radiation through
the same prepared heat network. For additional links supplied by another
physical owner, use
`prepare_with_directional_conduction_and_link_conductances` and continue to
provide their live values through `resolve_with_link_conductances`. Such links
must be disjoint from the oriented Fourier bindings. A source-integrated
isotropic layer cannot replace an explicitly directional link; it can still
serve other links in the same assembly.

This is a one-dimensional link/Galerkin reduction with a declared gradient
direction. It does not solve arbitrary three-dimensional anisotropic flow or
the transverse gradients of a generally oriented insulated rod. Those require
the spatial tensor-conduction owner. Source intervals, nonnegative principal
conductivities, unit directions and unchanged material endpoints remain explicit
admission conditions.

## Planar induction and finite heat

The prepared assembly also exposes `step_eddy_thermal_candidate` for the
original planar induction history. It uses the same thermal interfaces,
integrated material layers, directional conductivities and enclosure radiation
as ordinary heat/circuit stepping. All coefficients use the candidate cell
temperatures that drive the field's conductivity and magnetic storage law.
The original field solver supplies local Joule heat and parameter work; the
shared magnetic caloric owner supplies reversible heat and `U = psi + T*S`.

The caller supplies the original triangle-to-cell regions and source-backed
magnetic thermal laws. The prepared heat assembly remains the authority for
the material cards, caloric cells and links. Thermal contact resistance and
emissivity must be independent of unowned field coordinates. Graph resistors,
graph electrical contacts and circuit source bindings refuse because a spatial
field step does not evolve their current/work history.

The accepted report retains the complete `HeatNetworkSample` in
`report.coupled.source.material.network`, including per-contact, integrated
conduction and radiation results. Its top-level external-heat and physical
internal-energy accounts exclude the closed enclosure's roundoff net. Both
the native field/heat balance and the physical magnetic/caloric balance must
pass before either accepted history is returned.

`step_nonlinear_eddy_thermal_candidate` uses the same assembly for a nonlinear
field history. Each receiving region selects a `MagneticConstitutiveThermalLaw`:
the existing linear permeability model or an explicitly sourced separable
`H(B,T) = s(T) H_ref(B)` model. Linear conductors and nonlinear cores can share
one thermal solve. The positive affine scale and original reference curve stay
source-pinned; arbitrary two-coordinate magnetic laws are not inferred.
Backward-Euler magnetic convexity loss remains a numerical export and is
never deposited as material heat.

## Complete heat assemblies and material field storage

`contacts::heat::HeatAssembly` is the common interface for the prepared heat
network plus additional physical heat laws. It exposes the exact original
`PreparedHeatNetwork`, resolves one complete candidate state, and retains an
associated source sample that refers to its original `HeatNetworkSample`.
The ordinary prepared network implements this interface directly. Source
failures remain typed; cancellation and native solver/work refusals retain
their original error classes.

The airflow `SourceConvectionHeatNetwork` implements it for prescribed
reservoirs. Its `with_links` view adds source-resolved finite coolant using the
existing reciprocal native links. Prepare the underlying heat network with the
original initial additional-link conductances, then keep the same provider and
source laws across physical steps. The original complete/disjoint link and
exposed-area admission remains mandatory.

```rust,ignore
use fs_couple::circuit_enthalpy::field_storage::heat::field_storage_candidate;

let provider = convection.with_links(
    &finite_links, maximum_links,
    |binding, thermal_trial, property_temperature| {
        resolve_original_fluid(binding, thermal_trial, property_temperature)
    },
);
let proposal = field_storage_candidate(
    &provider, &cx, &thermal, target_s, &circuit, &ports,
    &field_bindings, controls, &energy_guess_j,
)?;
let (next_thermal, next_circuit, report) = proposal.into_parts();
let accepted_heat_and_fluid = &report.coupled.source.network;
```

This entry uses the original material field-storage driver. At each nonlinear
trial, the same finite-cell temperatures determine source-bound resistance,
C(T)/L(T), contact/layer conduction, radiation and fluid transport. The original
driver retains charge/flux, accounts for field-parameter work and field entropy,
and supplies reversible caloric heat separately from Joule loss. Changing only
fluid caloric or transport data can therefore change actual current, stored
field energy and subsequent temperature history.

Finite coolant exchange stays internal native conduction. Prescribed convection
is external heat. The original physical field-plus-caloric energy gate excludes
the retained closed-enclosure numerical net from external supply. No second
thermal solve, field-history advance or heat deposition is introduced.
`PreparedHeatNetwork::field_storage_candidate` delegates to this same entry,
so the existing complete-network path keeps its behavior and result defaults.

The field law still needs its original fully solid, fixed-volume material and
positive combined caloric/field capacity over its declared coordinate envelope.
Finite coolant cells keep their own caloric charts; they are not forced to use
the solid field cell's law. Pressure, composition and flow retain their declared
provider ownership. This does not add fluid momentum, advection, pressure work,
or new electrical/magnetic constitutive laws. Time refinement remains necessary.

## Reciprocal piezoelectric motion with the same heat assembly

`circuit_enthalpy::piezoelectric::heat::step_candidate` consumes the same
`HeatAssembly`, including the airflow provider and its finite-link view.
It advances the original rod, circuit and finite thermal history in one
candidate. The existing raw-coefficient piezoelectric API delegates to the
same private driver.

```rust,ignore
let proposal = fs_couple::circuit_enthalpy::piezoelectric::heat::step_candidate(
    &provider, &cx, &thermal, target_s, &piezoelectric_circuit,
    &electrical_ports, force_n, controls, &energy_guess_j,
)?;
let (next_thermal, next_model, report) = proposal.into_parts();
let accepted_heat_and_fluid = &report.source;
```

Source-bound load resistance, convection, contacts/layers and radiation are
resolved at the same proposed temperatures. The original reciprocal electrical
and mechanical solve supplies the actual current, force-displacement work and
Joule heat. A coolant material substitution can therefore change load
temperature, current, charge and continuing rod motion without changing that
solver or inserting a motion multiplier.

The entry checks accepted conductor identity and native heat destinations,
refuses unsupported electrical-contact work and unowned field-dependent heat
boundaries, and retains the complete accepted source sample. Closed radiation
remains internal under the original total-energy gate. Source failure, budget
exhaustion or cancellation returns no partially accepted thermal or mechanical/
electrical history.

This is the base small-signal rod model: ceramic coefficients, source context
and reference inertia remain independently held. The finite load and coolant
temperatures evolve. Temperature-dependent ceramic storage retains its
separate existing thermal entry; this adapter does not infer it from cooling
the load. Fixed heat geometry and the original first-order thermal coupling
require the same domain admission and timestep refinement.

## Temperature-dependent ceramic response with finite coolant

`circuit_enthalpy::piezoelectric::thermal::heat::step_candidate` accepts the
same `HeatAssembly` and an explicit `PiezoelectricThermalBinding`. Its
original ceramic driver resolves compliance, signed piezoelectric coupling
and permittivity at the same proposed temperature as the coolant and load.
The complete ceramic source, caloric receiver and rod reference volume must
agree. Changing the coolant transport or heat capacity can therefore change
the actual ceramic coefficients, charge and continuing mechanical response.

```rust,ignore
let proposal =
    fs_couple::circuit_enthalpy::piezoelectric::thermal::heat::step_candidate(
        &provider, &cx, &thermal, target_s, &piezoelectric_circuit,
        &electrical_ports, force_n, &ceramic_binding, controls, &energy_guess_j,
    )?;
let accepted_heat_and_fluid = &proposal.report().coupled.source.network;
```

The existing driver retains ceramic free energy, entropy and internal energy.
Reversible ceramic heat and native Joule heat enter their original caloric
destinations once. Finite coolant exchange remains internal; prescribed
reservoirs remain external, and the original physical energy gate excludes
closed-radiation roundoff from external supply. Typed provider failures remain
inside `PiezoelectricThermalError::Heat`; cancellation and work refusal retain
their native classifications. The concrete prepared-network method delegates
to this entry.

The finite-coolant ceramic regression compares continued endpoints with an
independent 2-by-2 constitutive inverse, midpoint circuit/motion equations and
caloric-plus-internal-energy balance. It changes viscosity, coolant heat
capacity and ceramic coefficient slopes only in source data, and checks
retained source identity and refused-candidate history. This uses the existing
fully solid, uniform, small-signal ceramic and positive combined capacity
envelope. Thermal strain, depoling and dielectric loss require their own
admitted laws; fixed heat geometry and first-order endpoint coupling remain.

## Thermoelectric conversion and field-actuator motion

The free functions
`contacts::heat::thermoelectric::step_thermoelectric_candidate` and
`step_thermoelectric_motion_candidate` consume the same `HeatAssembly`.
The original Kelvin-leg driver resolves integrated Seebeck voltage, signed
terminal Peltier heat and interior Thomson allocation at the actual trial
temperatures. Finite coolant now shares that solve with material Fourier
transport, load/contact resistance, thermal interfaces and radiation.

```rust,ignore
use fs_couple::circuit_enthalpy::contacts::heat::thermoelectric::{
    step_thermoelectric_candidate, step_thermoelectric_motion_candidate,
};
let proposal = step_thermoelectric_candidate(
    &provider, &cx, &thermal, target_s, &circuit, &ports,
    &kelvin_legs, controls, &energy_guess_j,
)?;
let moving = step_thermoelectric_motion_candidate(
    &provider, &cx, &thermal, target_s, &field_motion, &ports,
    force_n, &kelvin_legs, controls, &energy_guess_j,
)?;
```

The moving entry uses the original reciprocal electrical/field/mechanical
descriptor, including electrical back-reaction and native inertia/support
storage. Coolant material changes can therefore affect junction temperatures,
generated current and actuator motion through the same accepted endpoint.
Each entry returns one proposal from its supplied accepted history; the two
calls above illustrate alternative consumers, not sequential physical steps.

The original circuit deposits Joule loss once. Kelvin port work and reversible
junction/bulk heat remain paired internal transfers; finite coolant exchange
also remains internal. The physical energy gates retain actual external
electrical/mechanical work and prescribed heat, and exclude closed-radiation
roundoff. Full assembly samples remain in `report.coupled.source.network`.
`ThermoelectricNetworkError<E>` retains typed provider failures in `Heat(E)`
while its default preserves the concrete prepared-network API.

The finite-coolant Kelvin regression uses an affine Seebeck source and checks
independent four-cell heat balances, midpoint charge, reversible entropy and
total energy over continued steps. Only fluid viscosity or heat capacity
changes between cases. Existing Kelvin-leg geometry and single-phase limits,
source pins and refusal of unresolved contact temperature jumps still apply.
The field-motion entry retains its original held field/support material
snapshots; this connection does not add field entropy or moving heat geometry.

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

## Material-law feedback through the same heat assembly

`constitutive_graph::enthalpy::heat_assembly_candidate` connects the public
material-graph residual to the same complete `HeatAssembly`. Each proposed
thermal endpoint drives both the graph's temperature inputs and the assembly's
source-resolved transport. Original graph equations, exact material tangents,
nested Newton work bounds and native backward-Euler thermal storage retain
their existing owners.

```rust,ignore
let proposal =
    fs_couple::constitutive_graph::enthalpy::heat_assembly_candidate(
        &provider, &cx, &thermal, target_s, graph_input,
        graph_thermal_controls, &energy_guess_j,
    )?;
let accepted_heat_and_fluid = &proposal.report().source;
let (next_thermal, next_material_graph, report) = proposal.into_parts();
```

Finite coolant is reciprocal internal conduction. Prescribed reservoirs,
thermal contacts, integrated conductivity and radiation retain their own
source samples. Only explicitly bound graph outputs in watts become graph
heat; neither graph dissipation metadata nor a second circuit is silently
converted to a thermal source. Both prepared electrical bindings and nonempty
resolved resistor/contact work refuse. This preserves one owner for each
electrical equation and physical heat contribution.

The complete accepted provider sample remains in the graph thermal report.
Source failures and cancellation keep their original classes; other native
provider refusals retain their typed values in `GraphThermalError::Transport`.
The original closed-radiation correction excludes its numerical net from
external supply. Rejected trials publish neither thermal nor graph history.
The concrete `PreparedHeatNetwork::graph_step_candidate` delegates to this
same path.

A finite-coolant regression uses the original source-bound Ohmic graph node
under a solved voltage constraint. Independent current, two-body heat and
flat-plate convection equations check two continuing intervals under
viscosity, coolant heat-capacity and conductor-resistivity substitutions.
The graph law in that fixture is memoryless; the common driver retains its
existing general committed-history semantics. Fixed mass/geometry, admitted
source domains and first-order thermal coupling remain. Arbitrary graph
outputs still need their own complete work/storage account; this adapter
does not certify arbitrary material laws or add a chemistry model.

## Hysteretic cores and windings with finite material coolant

`ThermalHystereticCircuit::step_candidate_with_heat_assembly` uses the same
complete heat provider as other material circuit owners. Each trial starts
from the original accepted magnetic memory and resolves winding transport and
coolant at its proposed thermal endpoint. The native magnetic solve still owns
current, flux linkage, play history, stored magnetic energy and pinning loss.
Winding Joule heat and core pinning heat enter their original distinct caloric
cells once; finite coolant exchange remains internal.

A finite coolant link must already have a physical owner when the driver is
bound. Use `bind_with_link_conductances` with the initially resolved additional
link map. The original material network and that map must partition native
links completely and without overlap. These initial conductances admit and
resolve the initial state only: they are not cached as future cooling. Every
heat-assembly trial resolves the finite links again. The ordinary constructor
and bare step continue to refuse links without an owner.

The assembly must retain the exact bulk material network used at binding.
Full accepted fluid and heat samples remain in `report().source`, and
`HystereticThermalError<E>` preserves typed provider errors. Cancellation,
work budgets, immutable rejected candidates and the original closed-radiation
energy correction keep their existing semantics. The concrete
`step_candidate_with_heat_network` entry delegates to the same driver.

The finite-coolant regression checks two continued intervals against an
independent play-integral, winding/coolant balance and convection calculation.
Changing only coolant viscosity, coolant heat capacity or the source card's
pinning field changes the resulting current, magnetic memory and heating.
It also checks complete initial link ownership, fresh source identities and
queries, and transactional refusal after an accepted prefix.

The fitted magnetic coefficients remain constant over their explicitly
admitted temperature interval. Cooling changes winding R(T) and therefore the
magnetic trajectory; source queries refresh at the actual core temperature,
but no temperature-dependent magnetic free energy is inferred. Existing
fixed solid geometry, source coverage, monotone-linkage cuts and first-order
thermal coupling still apply.

## Biased magnetostrictive motion with finite material coolant

`ThermalMagnetostrictiveCircuit::step_candidate_with_heat_assembly` connects
the same complete heat provider to the original biased winding and reciprocal
rod-motion solve. Cooling changes winding temperature and its source-derived
resistance, which changes actual current, magnetic linkage and axial motion.
Every thermal trial starts from the same accepted linkage, extension and
momentum and resolves the finite coolant at that proposed endpoint.

Use `new_with_link_conductances` to admit links owned by the coolant provider.
The original material network and initial conductance map must cover native
links completely and without overlap. Initial values are admission data only;
a bare step still refuses links without their physical owner. The continued
assembly must retain the exact original bulk material network.

The existing biased work convention remains essential: the winding receives
`R(T) * (I_bias + delta_i_mid)^2 * dt`, including DC and cross terms, and
magnetic/mechanical storage change includes `I_bias * delta_lambda`.
These contributions are computed by the original driver and deposited once.
Finite coolant exchange is internal, closed radiation keeps its original
energy correction, and the full accepted provider sample is retained in
`report().source`. Typed provider failures, cancellation and work budgets
preserve both accepted histories.

The focused coolant fixture independently solves the reciprocal midpoint
current/motion equations and two-body heat balance over continuing intervals.
It changes only coolant viscosity, coolant heat capacity or the rod card's
piezomagnetic coefficient, checks actual motion and source identity, and
exercises initial ownership and later transactional refusal.

The rod's magnetic/mechanical coefficients remain at their independently held
preparation temperature and DC-biased reference. The finite winding temperature
does not infer core heating, thermal strain, hysteresis or a temperature-dependent
magnetostrictive law. Original amplitude bounds, stationary winding geometry
and first-order coupled thermal accuracy remain authoritative.

## Energy and scope

The standard circuit/enthalpy report preserves the actual source arithmetic,
including tiny net-enclosure roundoff. Radiation in this closed enclosure is
internal exchange, not external supply. The additional combined gate checks
`report.total_energy_residual_j + dt * sample.radiation_net_power_w()` against
the existing total-energy tolerance. True non-electrical external heat is
`report.non_electrical_heat_j - dt * sample.radiation_net_power_w()`. No surface
heat is clipped or redistributed to manufacture an exactly zero total.

This is fixed-geometry, fixed-topology Ohmic transport with uniform thermal
cells, supplied contact areas/heat splits, and a closed opaque gray-diffuse
radiation enclosure. Thermal and electrical contact resistances are independent
properties, not inferred from one another. View factors stay fixed. The ordinary
circuit entry holds capacitance/inductance fixed; the field-storage entry uses
their explicitly admitted source laws. Radiation may coexist with conduction or convection
on an exposed surface because these are distinct mechanisms. No participating
gas radiation, contact mechanics, spreading-temperature law, arc or flowing
melt is implied. Overall time accuracy remains that of the existing coupled
endpoint/midpoint method, generally first order. Energy balance is not an
experimental validation or a time-error certificate.
