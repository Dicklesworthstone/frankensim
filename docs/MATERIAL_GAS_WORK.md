# Material-driven closed gas work

`fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel` connects
source-backed species thermodynamics to conserved gas storage, heating and
reversible adiabatic compression/expansion. It owns one existing
`IsochoricGasParcel`; there is no parallel caloric law or gas-state registry.

```rust,ignore
use fs_material::thermochemical_gas::energy::{GasEnergySolve, entropy::GasEntropySolve,
    parcel::ClosedGasParcel};

// model is a ThermochemicalGasModel with caller-supplied NASA-9 species,
// composition and pure-species transport MaterialCards. point is a QueryPoint
// with explicit absolute temperature, pressure and any other source conditions.
let mut gas = ClosedGasParcel::new(model, point, 0.01)?;
let entropy = GasEntropySolve {
    minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
    entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128,
};
let energy = GasEnergySolve {
    minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
    energy_tolerance_j_kg: 1e-7, maximum_iterations: 128,
};
let compression = gas.set_volume_isentropic(0.0025, entropy, || false)?;
let heat = gas.add_heat(100.0, energy, || false)?;
let expansion = gas.set_volume_isentropic(0.01, entropy, || false)?;

// The updated common state is also the existing acoustic-medium interface.
let acoustic_medium = gas.state().state();
let pressure_pa = acoustic_medium.pressure;
let sound_speed_m_s = acoustic_medium.sound_speed;
let net_mechanical_work_j = compression.work_by_gas_j + expansion.work_by_gas_j;
```

The gas law, temperature-dependent heat capacity and entropy remain owned by
`fs-thermochem`. The material owner re-resolves density, pressure, viscosity,
conductivity, heat capacity, sound speed and acoustic impedance at each accepted
state. Changing composition or admissible species/transport source data changes
these responses without changing the heat or compression solver. All extra
query coordinates and source selections remain intact. A newly resolved state
is not automatically installed in an already compiled downstream solver; pass
it through that consumer's existing material/medium interface.

Compression conserves mass and specific entropy and uses `cv(T)/T` as its
thermal tangent. Heating conserves mass and volume and uses `cv(T)`. Both
inversions share bounded safeguarded Newton/bisection, explicit absolute
residuals and source-domain checks. Neither substitutes a frozen gamma or a
bracket-width-only convergence test. Failure and cancellation never publish
partial geometry or material state.

## Work and physical scope

Work is positive when delivered **by** the gas. For the admitted reversible
adiabatic path, `work_by_gas_j = -m*(u_final-u_initial)`. The result reports the
actual total entropy residual separately. The work value is a thermodynamic
endpoint calculation, not an independent mechanical-work measurement. The
regressions independently check it against integrated pressure-volume work.

This is a spatially uniform, closed, frozen-composition ideal-gas model on a
caller-admitted continuous stable branch. Heat and volume operations are
separate, not a time integrator for concurrent heating and moving boundaries.
There is no combustion/reaction, intake/exhaust, piston dynamics, shock,
friction, condensation, real-gas closure or experimental engine calibration.
A prescribed-heat idealized cycle must not be labeled simulated combustion.
An external mechanical consumer still owns exterior-pressure work, friction
and inertia; it must not count gas work twice.

Focused native regression command:

```sh
cargo test -p fs-material --lib thermochemical_gas::energy::
```

Tests include energy/entropy reconstruction agreement, the isentropic
`dp/drho = c^2` acoustic identity, nonlinear heat-capacity swaps, independent
pressure-volume work, a closed idealized cycle, unchanged-volume behavior,
source failures, iteration exhaustion and cancellation rollback. Synthetic
fixtures establish numerical behavior, not measured gas accuracy.
