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

## Same-card runnable consumer

`material_solid_self_heating` uses this public adapter with the original
material resolver and heat-capacity-curve constructor. Its insulated reference
specimen is the unit tetrahedron with vertices `(0,0,0)`, `(1,0,0)`, `(0,1,0)`,
`(0,0,1)` metres. Only the last node's z component is free. Its initial axial
coordinate is 0.001 m and velocity 0.1 m/s; all subsequent loads are zero.
Those are explicit fixture geometry, supports and initial conditions, not a
claim of a calibrated object. It writes physical motion, temperature, mechanical
storage, caloric-storage change, cumulative physical heat and energy residual
on one clock. It does not synthesize a sound from temperature or a loss preset.

```bash
cargo run -p fs-couple --example material_solid_self_heating -- \
  --demo 1 1 8 0.03125 /tmp/self-heating-baseline.csv
cargo run -p fs-couple --example material_solid_self_heating -- \
  --demo 3 1 8 0.03125 /tmp/self-heating-viscosity.csv
cargo run -p fs-couple --example material_solid_self_heating -- \
  --demo 1 2 8 0.03125 /tmp/self-heating-capacity.csv
```

After `--demo`, arguments are viscosity-time scale, heat-capacity scale, number
of steps, step duration in seconds, and a **new** CSV path. The manufactured
card contains density 2 kg/m3, Young's modulus 120 Pa, Poisson ratio zero,
proportional Kelvin time `0.01 * scale` seconds, and specific heat `1 * scale`
J/(kg K). These deliberately simple values are synthetic, not a real substance.
The example creates the entire card first; its consumer has no demo-only
coefficient path. The second command changes decay and heating. The third
keeps motion and dissipated joules identical while halving temperature rise.

For an existing canonical normalized material-card pack:

```bash
cargo run -p fs-couple --example material_solid_self_heating -- \
  --pack /path/to/material.matpack YOUR_KELVIN_TIME_PROPERTY \
  8 0.03125 /tmp/self-heating-source.csv
```

The selected model is explicitly **isotropic proportional Kelvin-Voigt**:
`C_viscous = time * C_elastic`, with no Maxwell branches in this example.
The time-property name is supplied, not guessed. Its value must be in seconds
and nonnegative; zero explicitly disables loss, while missing loss data refuses.
Density, Young's modulus, Poisson ratio and that time are resolved together at
300 K. The full native isotropic stiffness tensor and consistent mass are used;
there is no hand-built oscillator in the implementation. Reference mass comes
from the same source-density tetrahedral assembly that drives inertia and is
then used by the finite caloric cell.

The **same card** must supply `specific_heat_capacity` and `density` across the
fully-solid 300–350 K caloric interval. The existing caloric constructor owns
integration and chart admission. The mechanical source claims are pinned and
requeried at every actual thermal substep and event; changed values or unsupported
conditions refuse rather than quietly freezing temperature-dependent data.
Imported card identity and initial/final mechanical source identities appear
in the CSV header. No material-name selection or missing-data fallback occurs.

This remains a fixed-geometry, zero-thermal-expansion approximation, appropriate
only where the chosen frozen mechanical model is justified. The mathematical
Kelvin law is declared over 0..100000 rad/s by this bounded example; this is not
measured bandwidth validation. The existing propagation/source-domain limits
can reject an imported specimen or time step. Inputs are limited to 1..=4096
steps and packs of at most 16 MiB. The complete trajectory must succeed before
the output file is created; existing paths are never overwritten. File-writing
failures are not an atomic filesystem transaction.

`cargo test -p fs-couple --example material_solid_self_heating` defines three
additional regressions against independent scalar mass/stiffness/decay and
caloric equations, material-only substitutions, deterministic repeats, missing
or changed sources, explicit input limits and cancellation. Native compilation
and execution remain necessary to validate the authored Rust implementation.
