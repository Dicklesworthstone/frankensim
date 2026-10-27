# Physical tensor-solid loss heats finite material cells

`fs_couple::tensor_material::solid::thermal` connects the existing Tet4/Tet10
full-memory structural runtime to the existing enthalpy timeline. A changing
viscosity or Maxwell spectrum in the input data changes motion and physical
loss; the actual loss changes finite cell enthalpy and temperature. Changing
heat capacity changes the temperature rise, not a separately authored loss.

The mechanical owner already integrates the complete Kelvin and Maxwell
resistance channels along its exact held-force trajectory. This adapter routes
those independently integrated joules to one explicitly identified thermal
cell per original element. Recoverable Maxwell energy remains mechanical
storage. A work-minus-storage residual is never deposited as heat.

## Use the existing owners

Construct `TensorSolidBody` from `ViscoElements::project_basis`, then its normal
`history`. Create and activate the original `DepositionSimulation` with its
actual borrowed caloric curves and finite masses. The two accepted clocks must
match. Supply `ElementThermalCell` declarations in **connectivity order**, not
cell-ID order, and use the same exact caloric curve objects as the thermal
owner. IDs must be unique and every mechanical element needs a receiver.

```rust,ignore
use fs_couple::tensor_material::solid::thermal::{
    ElementThermalCell, SolidThermalConfig,
};
use fs_thermal::enthalpy::{CellId, StepConfig};

let receivers = [ElementThermalCell {
    cell: CellId(30),
    curve: &caloric_curve,
    temperature_interval_k: [300.0, 330.0],
}];
let candidate = history.advance_thermal_candidate(
    &thermal, &existing_observer, target_s, &held_nodal_forces_n,
    &receivers,
    SolidThermalConfig {
        thermal: StepConfig { max_dt_s: 0.001, energy_tolerance_j: 1e-10 },
        energy_tolerance_j: 1e-9,
    },
    &gate,
)?;
// Publish all three histories together only after inspecting the report.
let (history, thermal, observer, report) = candidate.into_parts();
```

Use `()` for the observer when the ordinary thermal links and boundaries are
sufficient. Otherwise its source, conductance, internal heat-transfer, timestep,
interval-source, event and step hooks are retained. They may use the existing
source-backed material transport owners. Mutable observer history must be
isolated by `Clone`, as required by the native observer contract.

Each actual thermal substep, including CFL reductions and event cuts, advances
mechanics exactly once for the same duration. Its integrated loss becomes an
average heat rate on that interval only; it never overwrites a boundary or
leaks into the next interval. Subsequent calls retain displacement, momentum,
all Maxwell coordinates and both accepted clocks. Held-force discontinuities
belong at call boundaries. Temperature, phase, caloric identity and reference
inventory are checked before stepping, after every accepted step and event.
Any late refusal discards the whole proposed mechanical/thermal/observer call.

## Energy boundary

The raw thermal audit counts generated heat in `external_heat_j`, because the
thermal kernel sees it as a source. The enclosing report removes only this
internal transfer to obtain `external_thermal_heat_j`. Its combined residual is

```text
change(mechanical kinetic + elastic + Maxwell storage)
+ change(caloric storage) - incoming caloric energy
- external thermal heat - mechanical port work.
```

Mechanical storage/loss/work, actual representable heat routing, and combined
storage each have an independent absolute residual check. All residuals are
reported, never repaired by changing a physical heat channel. Additional thermal
cells and their birth energy remain inside the native thermal accounting.

## Explicit approximation limits

This is constant-coefficient viscoelastic self-heating. Geometry, supports,
material stiffness, viscosity, relaxation times, reference mass and stress-free
strain stay fixed. Each declared temperature interval must justify freezing the
mechanical law; it is a caller admission, not an inferred empirical validity
range. The adapter does not manufacture temperature dependence, expansion,
phase-support changes or reciprocal thermoelastic heat. A partially liquid
receiver or a temperature outside that interval refuses. Spatial correspondence
and agreement between caloric and mechanical reference densities remain the
caller's specimen-binding responsibility; a cell ID alone cannot prove them.

Thermal transport retains the original explicit method, CFL, cell-uniformity
and spatial-network assumptions. Average interval heating preserves integrated
material heat, not the exact within-step temperature profile. Refine maximum
timestep and cell resolution independently. The original dense mechanical
capacity and projected-basis/strain-domain limits remain unchanged.

`cargo test -p fs-couple --test tensor_solid_thermal` defines six focused native
regressions: an independent damped-oscillator/finite-heat reference, data-only
viscosity/caloric substitutions, regional Maxwell memory with conduction and
all observer hooks, late temperature/phase/budget rollback, exact caloric/clock
binding, and Tet10 quadrature-to-element heat routing. These are synthetic
numerical fixtures, not measured-material validation.
