# Pump-selected circulation on the original material trajectory

`SourceConvectionHeatNetwork::step_pump_driven_relaxation_with_links` selects both
pressure and mass flow inside each finite-enthalpy trial using the existing
fixed-speed `PumpHeadCurve`. It shares the pressure-driven coolant, convection,
advection, hydraulic-heating and mechanical/electrical transaction; it introduces
no new material law, thermal state or time integrator.

## Original equipment and source materials

`PumpCoolantLoop` lives beside `PressureCoolantLoop` in
`fs_airflow::graph::thermal::medium::convection::deposition::source::heat::relaxation::pressure`.
It contains original sections and thermal links, pressure-law limits, a pump
curve and `PumpOperatingControls`. There is no separate prescribed pressure,
flow or convective speed. Section zero is the pump inlet: its resolved liquid
density converts head into pressure, and its kinematic viscosity must fit the
curve's explicit applicability range.

```rust,ignore
let report = network.step_pump_driven_relaxation_with_links(
    &cx, &mut owner, &mut thermal, &held_native_forces,
    electrical_input, thermal_controls, &energy_guesses_j, &gate,
    PumpCoolantLoop {
        sections: &sections, links: &wall_links, pump: &original_head_curve,
        pressure_budget, operating_point: pump_controls, maximum_links,
    },
    |section, trial, temperature| {
        original_liquid_laws[&section.cell].at_temperature(temperature)
    },
)?;
let accepted = report.transport_source.as_ref().unwrap();
let mass_flow = accepted.coolant.hydraulics.mass_flow_kg_s();
let supply = &accepted.pump.supply;
```

Initialize finite link coverage through the original
`*_with_link_conductances` APIs, as for pressure-driven coolant. Retain original
caloric materials, source selections, pressure and composition through all
trials. Full typed liquid objects accompany the final hydraulic result; generic
APIs cannot prove that independently supplied transport/caloric sources match.

`SeriesPressureFlow<ResolvedLiquidState>::solve_with_pump` also exposes the
fixed-state equipment/load calculation without requiring a thermal owner. It
reuses the SAME bounded matcher as parallel-liquid/common-solid cooling. Each
section supplies one checked snapshot; pressure candidates borrow snapshots,
not cloned material cards. Finalization moves original owners with the exact
matched evaluations, without replay. Every probe retains original support gates.

## Work and accepted history

Each thermal trial resolves bulk coolant properties, matches pump and losses,
derives speed from original local volume flow/area, evaluates declared bulk/film
convection, and transports that same mass through the original caloric residual.
Temperature therefore changes pump-selected flow and subsequent material response.

Local `dp_i * Q_i` is deposited once into the corresponding coolant cell. Pump
output is retained separately. Its discrepancy from independently summed loss
must fit the pressure criterion and interval energy tolerance.
`hydraulic_power_residual_w` reports this numerical discrepancy; it is not put
into a compensating heater. `external_thermal_heat_j` includes local hydraulic
loss as in the pressure-driven API. Do not add pump output again. Whole-run
energy audits can use reported pump work and retain accumulated matching error.

The payload comes from the accepted thermal trial. Failed brackets, source or
viscosity support, budget exhaustion or either cancellation authority publish
no proposed mechanical, circuit, material-memory, fluid, thermal or clock state.
There is no last-good operating-point fallback. The pump curve is held equipment
data, not an independently advanced dynamic state.

## Limits and computational cost

This is nominal quasi-steady closed circulation with common fixed density and
fixed inventories/geometry. A positive pressure bracket is required; no flow or
viscosity extrapolation occurs. Viscosity changes network loss, not an invented
pump-curve correction. Motor/shaft efficiency, speed dynamics, cavitation,
compression, variable-density expansion, changing immersion and moving geometry
are outside the model. Frozen mechanical exterior loading is not retuned from
coolant temperature. Refine the original first-order thermal split and upwind
spatial approximation independently.

Pump/network evaluations are bounded within each thermal source evaluation.
Each series solve is O(section count); source cards are not cloned per pressure
probe. Temporary evaluations still allocate. Structural modes and fluid
quadrature are retained; native solver costs remain. No measured hardware
speedup, allocation-free stepping or real-time qualification is claimed.

## Focused qualification

```sh
rch exec -- cargo test -p fs-airflow --test relaxation_convection finite::pressure::pump
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal::relaxation::temperature::heat
```

Three new native regression groups cover independent intersection/fixed-drive
parity, original sources, powered/unpowered continuation, local/global hydraulic,
RC and caloric balances, data-only substitutions, and refusal/rollback/retry.
Shared fixture temperature bounds now explicitly match their existing kelvin
queries; production source admission is unchanged. Native compilation, execution
and formatting remain UNEXECUTED: Rust/Cargo/rustfmt/RCH/DSR are absent here.
Correcting the previously reported shared caloric setup failure is not a native
pass claim for downstream tests.

Independent NumPy/SciPy references checked 18 static intersections and five
coupled cases at 256, 512, 1024 and 2048 Hz. Maximum relative intersection
difference was 2.46e-15; combined-energy discrepancy was 1.02e-10 J; individual
cell discrepancy was 1.12e-11 J. Refinement reduced errors against DOP853.
Doubling only viscosity changed synthetic endpoint flow from 0.193224 to
0.118059 kg/s and heated coolant from 303.952071 to 305.538367 K, with changed
charge and continuing displacement. These execute reference equations, NOT Rust,
and are not measured pump/material validation or a full-workspace pass.
