# Circulating coolant on the original coupled material timeline

`SourceConvectionHeatNetwork::step_advected_relaxation_with_links` joins the
existing prescribed mass-flow graph to the same finite coolant, Maxwell
material-memory, mechanical/fluid and optional circuit/field transaction.
A closed return loop transports energy internally. An inlet/outlet carries its
actual enthalpy supply through the existing finite cells. Neither is replaced
by a symmetric Fourier coefficient, an unlimited bath, or independently stepped
wall and coolant histories.

## One physical inventory and one endpoint residual

Prepare the original finite-coolant heat network and mechanical owner as in
`MATERIAL_FINITE_COOLANT_RELAXATION.md`. Define mass flows using those SAME
coolant `CellId` and caloric `MaterialId` values. The existing `Advected<()>`
owner admits prescribed fixed-cell mass balance and subtracts a canonical
common enthalpy datum before forming transport powers. Retain the flow graph
across calls; rebuilding it with different rates or inlet states is an explicit
caller intervention, not a fallback during a nonlinear solve.

```rust,ignore
use fs_thermal::enthalpy::advection::{Advected, MassFlow};

let flow = Advected::new((), vec![
    MassFlow { id: 1, a: coolant_supply, b: coolant_return, kg_s: mass_rate },
    MassFlow { id: 2, a: coolant_return, b: coolant_supply, kg_s: mass_rate },
], flow_budget)?;

let report = convection.step_advected_relaxation_with_links(
    &cx, &mut owner, &mut thermal, &held_native_forces,
    electrical_input, thermal_controls, &energy_guesses_j, &gate,
    &finite_convective_links, link_budget, &resolve_original_fluid, &flow,
)?;
let flow_at_endpoint = report.enthalpy_flow.as_ref().unwrap();
```

Use the native `MassBoundary::Inlet`/`Outlet` declarations for open circulation.
An inlet declares the original fluid material and supported specific enthalpy;
the outlet uses its current proposed cell enthalpy. Every flow cell must retain
balanced incoming/outgoing mass. The caloric storage must be enthalpy, not a
fixed-volume internal-energy chart masquerading as stream enthalpy. Existing
stationary-species admission remains with the original advection owner; this
entry does not evolve composition or replace material laws.

A bound structural material layer is NOT a flowing control volume. Such a path
refuses before motion because its retained material history belongs to that
fixed solid inventory. Advective links do not create surface contact area;
finite convection still uses the original native patches and material resolver.
Several walls and finite coolant cells can share the same circulation network.
Relative surface speed used by a convection correlation and mass rate used by
advection are explicit inputs: the caller must supply physically consistent
values. No hydraulic relationship is guessed from their numerical similarity.

At every proposed caloric endpoint, convection, donor enthalpy, original
conductivity/interface/radiation laws and electrical material properties see
the SAME state. The native circuit starts each trial from its accepted charge
and flux. Independently integrated Maxwell heat retains its original layers.
The existing Newton/enthalpy solve handles the complete residual; no flow
integrator, mass update, extra thermal state or endpoint heat correction exists.
Only after all existing source, solver, energy and cancellation checks succeed
do the material, circuit, mechanical/fluid and caloric histories publish.
The borrowed flow observer's audit is not advanced by discarded or accepted
endpoint queries; the returned report is the account for this operation.

## Energy and entropy accounts

`enthalpy_flow.endpoint_transfers` retains the actual reference-subtracted
internal powers; their sign need not match the mass-flow orientation.
They debit/credit the same native cell residual and are excluded from external
energy. `endpoint_sources_w` retains the inlet/outlet power view.
`external_enthalpy_j` is signed net stream supply and is ALREADY included in
`external_thermal_heat_j`. Do not add it again. Closed loops report zero stream
supply, including when they transport heat between cells at unequal temperatures.

The original caloric entropy account includes signed internal Q/T contributions
as UNCLASSIFIED at-cell transfer, not passive-conduction production.
`internal_at_cell_entropy_j_k` exposes that term and is ALREADY included in
`entropy.unclassified_source_at_cell_j_k`. This does not claim a stream-entropy,
mixing or second-law certificate. Reservoir convection, conductive links,
Joule/material dissipation, reversible field heat and closed radiation retain
their existing distinct meanings and original source payloads.

## Shared lower-level implementation

`FrameRelaxationThermalSystem::step_heat_network_with_advection` composes the
same flow with another complete heat-source callback. The callback resolves its
original prepared heat network and must not add this graph's fluxes itself.
Existing calls without advection retain the same numerical path with no
transfers; an empty flow graph is tested against that path.

At the native thermal layer, `Advected<()>::endpoint_transport` evaluates the
existing upwind graph without advancing it. `BackwardEulerStep::evaluate_with_transport`
and `evaluate_with_transport_and_entropy` consume `TransportStepCoefficients`
through the original storage/rate/acceptance implementation. The explicit and
implicit paths now use the SAME pairwise transfer routine. Local energy checks
prevent opposite errors from hiding in a zero global residual.

The flow evaluator uses bounded graph scratch and sparse power maps; internal
transfers have the original `max_links` ceiling. No new global flow matrix is
assembled. The original nonlinear solver and finite-difference budgets remain;
source queries and graph/result copies still cost work. Structural modes,
fluid quadrature and mechanical source projections are not rebuilt per caloric
trial. No native speedup, real-time or allocation-free claim is made.

## Boundaries and qualification

This is prescribed balanced flow through fixed control volumes, not a pump,
pressure, momentum or variable-mass solver. Composition, geometry and transport
source support stay explicit. The mechanical exterior-fluid loading is not
silently refreshed from coolant film temperatures. Fixed mechanical storage,
whole-layer caloric inventories and first-order material/thermal splitting
remain. Upwind spatial diffusion and time error need separate refinement; a
small algebraic energy residual is not a trajectory error certificate.

```sh
rch exec -- cargo test -p fs-thermal --lib enthalpy::backward_euler::transport
rch exec -- cargo test -p fs-airflow --test relaxation_convection finite::advection
```

Four native thermal and three coupled-circulation regression groups are authored.
They cover closed/open transport, large common caloric datums, empty-path parity,
local/global heat, circuit/material feedback, reversed flow orientation,
retained histories, invalid solid-flow or mass-balance declarations and late
refusal/cancellation/retry. Native compilation, tests and formatting remain
UNEXECUTED: Rust, Cargo, rustfmt, RCH and DSR are absent in this executor.

Independent NumPy/SciPy references ran directed 2/3/5-cell transport at four
time resolutions, inlet/outlet datum checks, and six complete mechanical/
fluid/Maxwell/RC/coolant scenarios at four resolutions. Maximum combined-energy
error was 2.02e-12 J in the linear transport reference and 7.26e-11 J in the
coupled reference. Refinement reduced temperature and mechanical differences
against separately integrated continuous equations. In one synthetic closed
loop, changing only mass flow from 0.2 to 0.8 kg/s changed the heated coolant
endpoint from 303.857564 K to 300.480677 K and altered charge and continuing
motion, with zero external stream supply. These are manufactured equation
checks, not execution of Rust, measured-coolant validation or full-workspace
qualification.
