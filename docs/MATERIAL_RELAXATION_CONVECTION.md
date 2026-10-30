# Source-fluid convection on the shared material trajectory

`fs_airflow::graph::thermal::medium::convection::deposition::source::heat::
SourceConvectionHeatNetwork::step_relaxation` connects the existing material
convection owner to `FrameRelaxationThermalSystem`. Original fluid transport
at the declared bulk/film temperature now changes actual layer heat, electrical
response and continuing material-memory/fluid-loaded mechanical motion.

This adds no convection correlation, material law, physical temperature store,
electrical integrator or mechanical time stepper. The same heat-network,
field/circuit, Maxwell and finite-enthalpy transaction remains authoritative.

```rust,ignore
use fs_airflow::graph::thermal::medium::convection::deposition::source::heat::
    SourceConvectionHeatNetwork;

let boundary = SourceConvectionHeatNetwork {
    heat: &prepared_material_heat_network,
    cells: &original_convective_cells,
    maximum_cells: cell_budget,
    resolve_fluid: |binding, trial, temperature| {
        // Select an already prepared original source by physical binding.
        // Retain pressure/composition/source pins; do not fit a new transport tuple.
        original_fluid_laws[binding.cell].at_temperature(temperature)
    },
};
let report = boundary.step_relaxation(
    &cx, &mut owner, &mut thermal, &held_generalized_forces,
    electrical_input, // Some for an electrically bound owner, even at zero drive.
    thermal_controls, &energy_guesses_j, &mechanical_gate,
)?;
let accepted_convection = report.transport_source.as_ref().unwrap();
```

The resolver receives the SAME proposed thermal state as conductivity,
interfaces, radiation, electrical R(T), and source C(T)/L(T) when bound. Actual
rho, heat capacity, conductivity and viscosity feed the existing dimensionless
correlation; sourced expansion is required by the existing natural-convection
path. Reservoir temperatures, relative speeds, characteristic geometry and
correlation applicability remain explicit inputs. No fluid-name branch or
material-independent heat-transfer coefficient is substituted.

The original thermal cells supply exposed area after contact coverage. A
nonzero authored ambient coefficient on a source-convecting cell refuses,
rather than adding duplicate cooling. Independent authored heating remains.
Convection is signed external reservoir heat, whereas enclosure radiation is
internal exchange. Solid relaxation/Joule losses and reversible field exchange
retain their own accounts. The generic caloric entropy report leaves the
callback's heat explicitly unclassified at the receiver; the retained source
payload contains its actual reservoir and transfer data. Q/T_cell alone is
not reservoir entropy or an irreversible-production claim.

`transport_source` is returned from the actual accepted thermal evaluation,
not a new query after publication. Source errors retain their original type.
Late source/solver/budget failure or either cancellation authority publishes
none of the electrical, material, mechanical, fluid, caloric or clock histories.
The initial heat assembly still must match original cell/link objects and the
accepted circuit's sources, geometry, interfaces and heat destinations.

## Shared composition rather than another solver

The lower-level `step_heat_network_with` accepts complete original heat
coefficients plus an additional typed source payload. `step_heat_network`
itself delegates through this same path with a unit payload. The source adapter
resolves the original prepared network once per thermal evaluation and adds
its physical exchange once. This preserves the crate dependency direction:
fs-airflow composes fs-couple, not the reverse. Original frame construction
types are reexported by `fs_couple::frame_vibration`; no alternate solid type
or additional package dependency is introduced.

The complete mechanical and thermal histories are reused. Structural modes,
fluid quadrature and source projections are not rebuilt for each caloric
trial. Material/correlation queries and the existing nonlinear solve still
have costs; no native throughput or real-time improvement is claimed.

## Physical boundary

These are PRESCRIBED reservoirs, not finite coolant inventories or solved flow.
This entry does not add finite-fluid thermal links, advection, fluid momentum,
phase change, moving contact area or thermal expansion. The existing first-order
rate/thermal split and whole-layer caloric model still require time/spatial
refinement. Correlation geometry and source support must be admitted separately.

The mechanical exterior-fluid type and convective-fluid type are independent
explicit bindings. The new heat path does NOT silently refresh frozen Stokes
added inertia/drag using a film temperature or a substituted coolant. Its
thermal feedback changes the admitted solid relaxation rates and continuing
motion; a live variable-density mechanical boundary needs its own state and
energy-transfer implementation.

## Qualification

```sh
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::modal::relaxation::temperature::heat
rch exec -- cargo test -p fs-airflow --test relaxation_convection
```

Two source-composition regressions and three concrete convection regression
groups are authored. They exercise accepted payload/parity, original coolant
identity and conditions, independent h and RC formulas, separate layer/global
energy, source-only substitutions, foreign history retention, empty-boundary
parity and typed late refusal/budget/cancellation with retry. Native Rust
compilation, runtime tests and formatting remain UNEXECUTED: Rust, Cargo,
rustfmt, RCH and DSR are absent in the editing executor.

Independent NumPy/SciPy reference equations passed six powered/unpowered,
wet/dry, warming/cooling cases at four timestep resolutions and a stationary
heated control. Maximum combined-energy discrepancy was 7.05e-11 J. Refinement
reduced temperature and trajectory errors against continuous DOP853. Changing
only coolant conductivity in one synthetic warm-reservoir case changed the
coating endpoint from 301.888088 K to 304.066077 K and changed both charge and
continuing displacement. A fully relaxed stationary specimen stayed motionless.
These manufactured equations do NOT execute Rust, validate a real coolant or
specimen, or qualify the geometric correlation or full workspace.
