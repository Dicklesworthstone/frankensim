# Material-driven electrical heat, phase storage and gas motion

`GasPiston::powered_thermal_candidate` solves a state-dependent body power
alongside the existing finite-body heat exchange and gas/momentum equations.
It is a public library capability, not a prescribed-watts example:

- `W = dt * P(gas_mid, body_mid)`;
- `Q = dt * G(gas_mid, body_mid) * (T_body_mid - T_gas_mid)`;
- the body stores its prescribed heat plus W minus Q;
- the gas receives its prescribed heat plus Q and exchanges work with the piston.

Every proposal starts from the same accepted state. Both midpoints come from
conserved energy coordinates, not an average of endpoint temperatures. The
existing bounded secant/bisection owner solves source and exchange feedback;
the existing piston and gas owners retain momentum, pressure-work and caloric
inversion. The final audit uses actual evaluated source power, not merely the
trial energy supplied to the inner solve. Source, exchange and total-energy
residuals have explicit absolute SI tolerances and independent iteration caps.

## One source binding for conduction and Ohmic heating

`thermal::power::joule::MaterialJouleHeating` calls the existing
`fs_material::conductor` resolver, ideal-current/voltage law and Joule-power
operation. It does not select a metal, copy a resistivity coefficient or own a
second resistance law. The original resistivity and thermal properties are
queried through the same `ThermalMaterialPoint`, bound to the same card as the
body's enthalpy chart. Temperature and explicitly declared liquid mass fraction
change; all other source coordinates and quantity descriptors remain intact.

```rust,ignore
use fs_couple::cavity_phs::piston::thermal::{
    source::ThermalMaterialPoint,
    power::joule::{MaterialJouleHeating, JouleHeatingError},
};
use fs_material::conductor::OhmicDrive;

// Existing card, context, caloric body input, keys, selections and controls.
let heater = MaterialJouleHeating {
    material: ThermalMaterialPoint {
        material: &card, context: &context, temperature_axis: "T",
        liquid_fraction_axis: Some("melt"), // must already exist in context
    },
    selection: electrical_selection,
    length_m: 0.01, area_m2: 0.000025,
    drive: OhmicDrive::Voltage(50.0),
};
// Shares heater's exact source conditions; heat-transfer area is NOT wire area.
let layers = heater.conduction_layers(&conductivity_key, thermal_selection,
    0.02, 0.001, 0.02);
let proposal = piston.powered_thermal_candidate(
    input, body, piston_controls, exchange_controls, power_controls,
    |gas, phase| layers.resolve(gas, phase).map_err(JouleHeatingError::Thermal),
    |_, phase| heater.resolve(phase, body.volume_m3),
    || gate.is_requested(),
)?;
let (coupled, electrical, power_balance) = proposal.into_parts();
let (next_piston, next_body, heat_balance) = coupled.into_parts();
// Publish both endpoints together. Electrical heat and exchanged heat are
// ALREADY present in these states; never add either a second time.
```

The electrical geometry must represent the whole uniform thermal body
(`length*area == body.volume`, to roundoff). Resistivity uses the existing
`electrical_resistivity` scalar SI schema. Thermal conductivity retains its
original dimensional or static semantic key. Missing data, mismatched source
identities, unsupported coordinates and dynamic/frequency-only laws refuse,
even with zero drive. Explicit zero drive switches off electrical input without
resetting thermal, phase, gas or mechanical history. Source records retain the
actual conductor, signed current/voltage and accepted body midpoint.

Increasing resistance increases power under ideal current control and reduces
power under ideal voltage control. During a latent plateau, a source-declared
phase-fraction resistivity can change both power and melt fraction without
inventing a temperature increase. No implicit alloy or mixture rule is inferred.

## Scope and focused checks

This is a fixed-geometry, spatially uniform Ohmic body, a closed frozen-composition
gas and the existing translating wall. It is not a spatial current solver,
contact resistance, AC circuit, magnetic coupling, electrochemistry, combustion,
flowing melt or thermally deforming piston. A finite electrical supply requires
its own coupled energy-store law; the drives here are ideal prescribed sources.
Source brackets and timestep refinement remain caller obligations. Convergence
refusal is not a theorem of runaway, and no unconditional stability, monotonicity
or nonnegative entropy-production guarantee is made.

Eight new Rust tests cover original source receipts, current/voltage signs,
resistivity-only changes through full trajectories, latent-phase feedback,
shutoff/cooling, fixed-power parity, an independent analytic feedback solution,
invalid source/geometry/static support and cancellation/refusal atomicity.
The existing instrument-material workflow already selects this module:

```sh
cargo test -p fs-couple --lib cavity_phs::piston:: -- --test-threads=1
```

Synthetic equation checks are not measured-material qualification or a native
Rust test pass. Local Rust execution was unavailable during authoring.
