# Workload transients with coupled cooling

The experimental `cooling-network` command can advance a heat-storing solid
through a piecewise-constant workload and positive fan-speed schedule. It uses
the existing tetrahedral conduction/contact operator and mixed-air network.
It is not a native `.fsim` project, transient CFD, or a ledger-backed run.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/transient-contact-pulse.json
cargo test -p fs-conduction --lib transient::backward_euler
cargo test -p fs-cli --bin frankensim network_command::transient
cargo test -p fs-cli --test transient_cooling_network
```

## Request

The optional top-level `transient` object selects time integration rather than
the original steady solve. It requires `objective.gradient=false` and excludes
both `design` and `fan_speed_design`: steady derivatives and steady design
criteria are not silently reused for a transient trajectory.

```json
"transient": {
  "initial_temperature_k": 300,
  "volumetric_heat_capacity_j_m3_k": 2000000,
  "max_step_s": 2,
  "max_steps": 1000,
  "temperature_limit_k": 305,
  "intervals": [
    {"duration_s": 30, "power_scale": 1, "fan_speed_ratio": 1},
    {"duration_s": 120, "power_scale": 0, "fan_speed_ratio": 1.5}
  ]
}
```

Choose exactly one initial-temperature representation: the uniform scalar
`initial_temperature_k`, or `initial_temperatures_k`, an array with one positive
kelvin value per solid vertex. Likewise choose uniform
`volumetric_heat_capacity_j_m3_k` or `element_heat_capacities_j_m3_k`, one positive
J/(m³ K) value per tetrahedron in input order. Capacities are explicit caller
declarations, not inferred from conductivity or material names.

`intervals` is nonempty, starts at zero, and contains consecutive positive
durations. Each interval selects exactly one of a nonnegative `power_scale`
on the complete original source field, or `component_powers_w`, an object giving
absolute nonnegative watts for every declared component. These modes never
combine or inherit missing values from an earlier interval. A fan-driven request must supply a positive admitted speed in every
interval. A pressure-driven request must omit `fan_speed_ratio`. A zero-speed
fan and a zero-flow exchanger remain unsupported: this does not invent a
natural-convection fallback for a stopped fan.

Each interval is divided into `ceil(duration_s / max_step_s)` equal steps,
subject to floating-point endpoint representation. Thus no step crosses a load
or speed discontinuity. The complete step count is admitted before physics;
`max_steps` is capped at 10,000. The existing wall budget covers the entire
numerical invocation, not each time step separately. Interval coefficients,
air properties and imposed supplies are held constant. Fan/graph flow and all
flow-derived convection coefficients are recalculated for every interval.

## Discretization and coupling

The solid uses the consistent P1 capacitance matrix, with optional distinct
capacity per element, and backward Euler:

`(C + dt K) (T_new - T_old) = dt (b - K T_old)`.

The linear solve acts on the temperature change, avoiding an unnecessarily
large absolute-temperature right-hand side. Contact resistance and heterogeneous
conductivity remain in the spatial operator. Every Robin-reference coupling
iteration uses the **same previous accepted field**, not the previous coupling
trial. A converged endpoint becomes history only after solid, air and energy
checks pass. A failed step publishes no partially advanced state.

Air is quasi-steady at the endpoint, including downstream heating, mixing,
reverse-flow ordering and explicit bypasses. There is no air heat capacity,
transit delay, hydraulic inertia, fan acceleration or continuously ramped input.
The approximation needs a time scale on which treating the air as quasi-steady
is meaningful; the command does not certify that assumption for the user.

For each time step, the actual published temperatures are used to check
`stored_change = dt * (source_power - external_air_heat_gain)`.
Internal contact heat is not counted twice. Whole-window input, exhaust gain,
storage and the raw residual are reported as well. The absolute per-step joule
gate is `tolerances.heat_w * dt`; no time-integration error bound is inferred.

## Outputs and temperature limits

The original result field is the **final accepted transient field**, with its
matching fan speed, flows, coefficients and contact results. `solid_inputs`
contains the base declarations; history records either the applied `power_scale`
or the complete absolute `component_powers_w` map. The other field is `null`. The `transient` result includes time, total steps and solid solves,
energy totals, every accepted endpoint's objective and heat/storage values,
and the peak sampled objective with its sample time. The initial condition is
included in the peak calculation. Full nodal history is not retained.

`temperature_limit_k` is optional and monitors the selected mean or discrete
maximum objective. `first_sampled_violation_s` is the first recorded state above
the limit, **not an interpolated or certified crossing time**. The sampled peak
is not an upper bound between steps or on the spatial continuum. Halving the
step is an explicit resolution study, not an automatic certification. No
transient gradient, adaptive time stepping or durable checkpoint is supplied.
A final nodal field may be declared as another request's initial state; its
local timeline starts at zero again.

## Pulse example

The example applies 20 W to a localized component for 30 seconds, then removes
the input for 120 seconds and raises the fan ratio from 1 to 1.5. Two materials,
separately owned interface traces, finite contact resistance and two different
heat capacities are retained. With a 2-second maximum step, an independent
NumPy P1/contact/air reference gives a sampled hotspot of 306.3431585 K at
30 seconds and a final maximum of 301.4672336 K. A 305 K limit is first exceeded
at a recorded endpoint at 22 seconds: considering only the final field would
miss that earlier excursion. The 600 J input divides into approximately
536.2871881 J stored and 63.7128119 J passed to the external air.

FrankenSim itself (measured 2026-09-25, release build, after 05db922bf made the transient capacitance row-sum lumped; the reference above uses the consistent P1 mass, so the two differ on this 12-tet mesh) gives a sampled hotspot of 304.1166334 K at
30 seconds and a final maximum of 301.7305465 K. The example therefore now
declares a 303.25 K limit, keeping the limit's role. That limit is first
exceeded at 24 seconds, and the 600 J input splits into 533.2318718 J stored
and 66.7681282 J to air. The Rust tests pinning these values were executed
natively on 2026-09-25.

These are independent mathematical reference values, not retained executions
of the Rust implementation. Compilation, formatting and Rust tests have not
been run in the authoring environment. Backward Euler is first-order in time;
the supplied validation also compares refinements against the exact matrix
exponential of the same semi-discrete linear model, not the continuum PDE.

## Independent component workloads

A workload can move between fixed component footprints without replacing the
mesh or scaling every component together. Declare the footprints once in
`solid.component_power`, including a zero nominal power for initially idle parts.
Then each interval supplies **every component**, including explicit zeros:

```json
{"duration_s":30,"component_powers_w":{"cpu":20,"gpu":0},"fan_speed_ratio":1}
{"duration_s":30,"component_powers_w":{"cpu":0,"gpu":20},"fan_speed_ratio":1.5}
```

These are alternative interval objects, not extra top-level declarations.
Missing, extra and negative powers refuse; unknown names are checked across the
whole schedule before any numerical work. Supplying `power_scale` together with
`component_powers_w` also refuses. Named workloads require component footprints,
so a uniform `source_w_m3` cannot be silently split into invented components.

The existing `PowerMap` projects the new watts onto the same nodal P1 supports
once per interval. It does not divide by nominal wattages or try to recover
individual sources from an already combined field. Overlapping footprints still
superpose, and a component with zero nominal watts can turn on normally. The
projected and subsequently assembled powers are both checked. Only the current
interval's source field is stored, not one full nodal vector per interval.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/transient-component-workloads.json
```

This example moves a 20 W pulse from the spreader-side CPU footprint to the
substrate-side GPU footprint at 30 seconds, then turns both off at 60 seconds.
The named parts are illustrative labels, not a resolved semiconductor model.
Independent NumPy calculations place the sampled peak at approximately
319.642150 K at 60 seconds, with 1200 J total input. The hottest vertex moves
from 4 to 13. From the same cold initial condition and fan speed, a single 20 W,
30-second pulse produces about 306.343159 K at the CPU footprint versus
318.415974 K at the GPU footprint: equal total power is not equal hotspot risk.
These are independent mathematical references, not executed Rust measurements.
