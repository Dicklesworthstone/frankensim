# Finite coolant on the shared material trajectory

`SourceConvectionHeatNetwork::step_relaxation_with_links` connects the existing
finite-fluid convection owner to `FrameRelaxationThermalSystem`. A coolant has
its own ORIGINAL caloric cell, mass and material law, not an unlimited prescribed
temperature. Several structural surfaces can exchange with the same coolant.
Its changing temperature feeds back into material-resolved heat transfer, solid
relaxation, electrical transport and continuing mechanical/fluid-memory motion.

This adds no correlation, material evaluator, circuit solver or thermal/mechanical
integrator. The original thermal links apply equal and opposite energy, while
all accepted histories share the existing transaction boundary. Heat-capacity,
conductivity, viscosity or inventory changes are explicit data substitutions.

## Prepare one inventory and each physical heat path once

Define existing native `ThermalLink` objects for wall/fluid patches. Their areas
and per-side coverage remain authoritative; their authored conductance must be
zero. Exclude them from the material network's Fourier bindings. Keep the same
coolant card/context behind its caloric law and its prepared transport resolver.
The generic transport trait cannot independently prove that correspondence.

Resolve the initial links using `evaluate_source_link_convection`, then use the
returned conductance map for cold admission:

```rust,ignore
let initial = evaluate_source_link_convection(
    &thermal, &fluid_links, link_budget, &resolve_links, || gate.is_requested(),
)?;
let mut owner = FrameRelaxationThermalSystem::new_shared_with_link_conductances(
    material_bindings_and_layer_cells, mechanical, &thermal, &thermal_network,
    initial.conductances_w_k(), heat_budget, &gate,
)?;
let heat = heat_network.prepare_with_link_conductances(
    &cx, &thermal, initial.conductances_w_k(),
)?;
let network = SourceConvectionHeatNetwork {
    heat: &heat, cells: &reservoir_boundaries, maximum_cells: reservoir_budget,
    resolve_fluid: &resolve_reservoir,
};
let report = network.step_relaxation_with_links(
    &cx, &mut owner, &mut thermal, &held_native_forces,
    electrical_input, thermal_controls, &cell_energy_guesses_j, &gate,
    &fluid_links, link_budget, &resolve_links,
)?;
```

For electrical assemblies, cold-bind through
`bind_material_circuit_with_link_conductances` or
`bind_material_contact_circuit_with_link_conductances`. They reuse the original
bulk/interface admission and heat destinations. Their `fields` argument is empty
for authored fixed C/L, or contains complete original `FieldThermalBinding`
objects for material-dependent field storage. The ordinary initial binding
still refuses missing Fourier coverage instead of discarding a physical link.

Initial conductances are NOT cached boundary laws. Every proposed endpoint uses
both actual finite wall/fluid temperatures and the explicitly selected bulk or
film property convention. The same prepared fluid laws supply rho/cp/k/mu and,
when required, thermal expansion. No material name selects a correlation.
Ordinary stepping without the external link owner refuses incomplete coverage;
it cannot silently turn a finite coolant path into insulation or a fixed bath.

## Continuing state and physical accounts

`report.transport_source` contains `FiniteRelaxationConvection`: original
prescribed-reservoir results in `reservoirs`, and actual finite-interface source
results in `finite`. These are the final accepted evaluations, not a fresh
query after publication. Material objects, physical context, requested property
temperature, signed heat and native link IDs remain inspectable.

Finite convection enters the native conductance map, NOT a duplicate pair of
heater powers. It is internal energy exchange and contributes once to the
existing conductive entropy account. Only real prescribed reservoirs contribute
external convective heat. Empty reservoir bindings admit a closed finite thermal
network. Several walls may heat/cool one fluid cell through distinct patches.

Joule/interface loss, mechanical relaxation, reversible field heat and exterior
mechanical losses retain separate accounts. No Stokes loss is silently deposited
into a new bath. Failed source, coverage, budget, electrical or thermal admission,
and cancellation, publish no solid/coolant temperatures, charge/flux, material
histories, structural/fluid momenta or clock. Deterministic retry starts at the
same accepted state. Prepared source pins and native resource bounds remain.

## Model and performance boundaries

Finite caloric inventory does not imply solved coolant flow. Relative velocity,
pressure, composition and geometry are prescribed. This adds no advective mass
transport, fluid momentum, pressure work, phase-front motion or moving interfaces.
The coolant and mechanical exterior-fluid bindings remain independently explicit;
film-property updates do not silently rebind frozen Stokes inertia or drag.
Solid storage moduli/reference density remain fixed in this rate-only path.
The original first-order mechanical/thermal split and thermal lumping still need
time/space refinement. Correlation support and geometry validity remain explicit.

Structural modes, fluid quadrature and mechanical source projections are reused.
The new adapter delegates one complete source evaluation per thermal trial; it
adds no independent iteration or accepted thermal store. Source queries, native
rate rebuilds and the existing nonlinear solve still cost work. No native speedup
or real-time performance claim is made.

## Focused verification

```sh
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal::relaxation::temperature::heat
rch exec -- cargo test -p fs-airflow --test relaxation_convection
```

Four new native regression groups cover extra-link cold admission, omitted-link
refusal, empty-map parity, two walls sharing a finite coolant, powered/unpowered
and mixed-reservoir trajectories, local/global energy, conductive entropy,
reversed link orientation, data-only transport/capacity/inventory changes,
retained foreign material history, and typed late refusal/cancellation/retry.
Native compilation, execution and formatting are UNEXECUTED in this editor:
Rust, Cargo, rustfmt, RCH and DSR are unavailable. Source inspection and reference
arithmetic do not establish native build or test success.

Independent NumPy/SciPy equations passed seven wet/dry, powered/unpowered,
warming/cooling and mixed-reservoir cases at four timestep resolutions. Maximum
combined-energy discrepancy was 9.32e-11 J; maximum local-cell discrepancy was
1.19e-11 J. Refinement reduced error against separately integrated continuous
DOP853 equations. In the manufactured baseline the finite fluid cooled from
310 K to 307.812747 K while the coating warmed from 300 K to 300.826970 K.
Doubling only coolant conductivity gave 306.867668 K and 301.182612 K,
respectively, changing both charge and continuing motion. These synthetic
references do NOT run Rust, validate a real specimen or establish correlation
accuracy, a continuous error bound, full-workspace readiness or performance.
