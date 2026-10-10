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

## Fixed-density total enthalpy and latent heat

`transient.enthalpy` selects the spatial total-enthalpy owner instead of the
temperature/capacity formulation above:

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/enthalpy-phase-pulse.json
cargo test -p fs-cli --test cooling_enthalpy
```

The complete [phase-pulse request](enthalpy-phase-pulse.json) declares one
equilibrium phase chart for the solid. Each chart knot gives
`specific_enthalpy_j_kg`, positive `temperature_k` and
`liquid_mass_fraction` in `[0,1]`. The existing material owner validates the
2–4096 knots and evaluates the piecewise-linear enthalpy relation. Equal
temperatures over a positive enthalpy interval retain an exact latent plateau;
no apparent heat capacity or extra latent-energy source is added.

Inside `enthalpy`, provide positive `reference_density_kg_m3`, a nonempty
`source`, a 64-hex `material_card_identity`, and exactly one initial state:
`initial_specific_enthalpy_j_kg` or nodal
`initial_specific_enthalpies_j_kg`. The identity is a caller declaration,
not proof that a measured material card was resolved. The example says so
explicitly. Omit the temperature initializers, volumetric-capacity fields and
`transient.nonlinear`; mixing storage formulations refuses. Conductivity
declarations and supported finite thermal contacts remain independent inputs.

To declare an ordinary sensible solid without inventing a melting endpoint,
add `"phase": "solid"`; every knot must have liquid fraction zero and strictly
increasing temperature. `"phase": "liquid"` likewise requires fraction one.
The default `"phase": "solid-liquid"` retains the mixed-phase constructor,
including its fully solid/fully liquid endpoint requirements and exact plateaus.

For multiple storage materials, replace the uniform chart fields inside
`enthalpy` with a `materials` array and an `element_materials` array. Each of
the 1–256 material records supplies a unique `name`, `material_card_identity`,
`source`, `reference_density_kg_m3`, `knots`, and optional `phase` as above.
`element_materials` contains exactly one of those names for each solid
tetrahedron, in the same order. Initial nodal enthalpies and Newton controls
remain beside these arrays. See the complete
[two-material contact request](enthalpy-contact-materials.json).

Storage names never select conductivity or contact cards. Their assignments
are independent of `solid.materials` and `solid.element_materials`, and no
association is inferred from matching names. A single vertex cannot carry two
different storage identities, even when their tabulated values agree. Use
separate interface vertices and an explicit `solid.contacts` declaration for
different enthalpy charts; undeclared mixtures refuse. Shared vertices within
one storage material remain valid.

The required `enthalpy.newton` object declares `max_iterations`,
`residual_rtol`, `residual_atol_j`, `linear_restart`, `max_linear_cycles`,
`armijo_c`, `shrink` and `max_backtracks`. `linear_restart` is capped at the
smaller of the solid vertex count and 256 to bound Krylov workspace. The worst-case product
`max_iterations * linear_restart * max_linear_cycles` must fit
`budgets.linear_iterations`. Existing interval, step-count and wall-time
budgets apply. The endpoint residual is
`M_ref (h_new-h_old) + dt [A(T(h_new)) T(h_new)-b]` in joules, with lumped
reference masses `rho_ref * V / 4` per tetrahedron vertex. Its target is
`max(residual_atol_j, residual_rtol * norm(R(h_old)))`; it also must pass the
independent physical energy gate. All air-reference and radiation trials
reuse the same previously accepted enthalpy. Only a complete accepted endpoint
advances time or enters the accumulated energy account.

Ambient radiation uses the existing optional `radiation.surfaces` declaration.
It can remove heat or add heat from a hotter reservoir. The air network receives
only convection. Radiation has its own watt and joule fields, so the window
balance is `stored_change = input - air_gain - radiative_loss`.

The final result adds `solid_specific_enthalpies_j_kg` and
`solid_liquid_mass_fractions` beside `solid_temperatures_k`.
`transient.scheme` becomes `backward-euler-total-enthalpy`;
`transient.enthalpy` retains the declared chart identity/source/reference
density, initial and final total reference enthalpy, and solver controls/work.
For assigned materials it reports the material table with chart identities,
`element_materials` and the derived `vertex_materials` instead of the uniform
identity/density fields. Summaries use each vertex's own constitutive chart
and reference mass, including when only some materials melt.
Each history row also reports the minimum and maximum specific enthalpy and
reference-mass-weighted mean liquid fraction. Full nodal enthalpy history is
not included in the forward JSON. A new fixed schedule can start from the returned final nodal
enthalpies; temperature alone cannot reconstruct a state on a latent plateau.

The synthetic example stays at 350 K while its enthalpy and liquid fraction
change. A 0.4-second pulse supplies 400 J, followed by 0.4 seconds without
internal heating. Independent tetrahedral/heat-exchanger algebra gives
49.5595805903 W of convection throughout. With a 300 K black reservoir,
radiation removes 741.247808319 W and the final mean liquid fraction is about
0.43020623. Changing only that reservoir to 400 K gives a radiative loss of
−1137.02718109 W and a final mean liquid fraction of about 0.88099222; air
convection stays the same. The binary tests compare individual nodal enthalpies,
every accepted endpoint summary, air outlet temperature and accumulated energy
with these independent formulas, and compare a split/restarted trajectory with
the uninterrupted run. These are numerical references, not measured material
or hardware validation.

The two-material example has matching contact traces between mirrored unit
tetrahedra. Their reference densities are 10 and 20 kg/m³ and their latent
temperatures are 350 and 330 K. A 0.5 m² contact with resistance
0.1 m² K/W transfers 100 W from the first material to the second. This transfer
changes individual nodal enthalpies and cancels from the combined energy
balance. Two independent unit-capacity-rate air streams remove a total
78.0845054746 W. Radiation to 300 K reservoirs removes 902.815209140 W;
changing the reservoirs to 400 K adds 2059.88235101 W. The 0.2-second pulse
supplies 400 J, followed by 0.2 seconds of cooling. All nodes remain inside
their declared plateaus, so the binary test can check every nodal enthalpy
against independent surface-area/contact/reference-mass arithmetic, along
with mass-weighted phase summaries and split/restart parity. Additional
tests replace one chart by an ordinary solid and compare separate uniform
solid/liquid diffusion cases with an independent four-node linear solve.

This opt-in CLI mode supports fixed schedules and fixed reference density with
uniform or explicitly assigned equilibrium charts, including single-phase
charts. Geometry, mass and energetic internal variables are frozen. Adaptive
and repeated schedules, steady design, time/mesh studies, recirculation and
enclosure radiation explicitly refuse. There is no fluid storage, phase
advection, melting-driven motion or certified inter-step peak.

### Enthalpy history gradients and workload/fan sizing

Add `"adjoint": {"qoi": "sampled-peak", "max_checkpoint_bytes": 1048576}`
inside `transient` to differentiate the earliest sampled maximum. Set `qoi`
to `"final"` for the final objective instead. The dedicated enthalpy history
adjoint includes the coupled air references and the physical ambient-radiation
feedback. It carries the enthalpy derivative directly between steps, including
through latent-plateau interiors where temperature alone cannot carry the
stored-state sensitivity. It supports the same explicit single-phase and
heterogeneous chart assignments as the forward solve.

The `transient.adjoint` result includes:

- `dtemperature_dinitial_specific_enthalpies_k_kg_j`: one derivative per initial
  nodal enthalpy, in K kg/J.
- `dtemperature_duniform_initial_specific_enthalpy_k_kg_j`: the derivative for
  an equal additive change to all initial nodal enthalpies, in K kg/J.
- `intervals[].dtemperature_dpower_multiplier_k`: the derivative of a multiplier
  on that interval's actual declared workload. A zero-workload interval has
  zero derivative for this multiplicative control.
- `dtemperature_dinlet_temperatures`: derivatives of inlet temperatures applied
  throughout the trajectory, in the existing graph-node order.
- `intervals[].dtemperature_dlog_fan_speed_ratio_k`: the total speed derivative
  for the admitted single affinity-scaled fan bank. It includes changing air
  capacity, complete solid/air/radiation feedback, and supported Reynolds
  dependence of a flow-derived convection law. A declared constant HTC has
  zero convection contribution. This field is `null` without a fan or when a
  convection card lacks an admitted smooth Reynolds derivative.

The tape bounds retained h, temperature and air-reference checkpoints and
control accumulators with `max_checkpoint_bytes`, then reconstructs accepted
steps with exact h/temperature replay and renewed physical residual/energy
checks. This cap does not include all transient solver workspace. Derivative
work shares the original wall deadline; a failed replay or derivative solve
refuses the result. Chart slope corners and validity endpoints refuse classical
endpoint derivatives; ties retain the selected branch without claiming a
unique derivative. Geometry, chart data, reference density, conductivity,
contact resistance, fan curve, quadratic loss coefficients, convection law
data and radiation controls stay
fixed. The optional `component_power` and `contact_resistance` adjoint requests
are unsupported in this mode.

The complete [enthalpy power-sizing request](enthalpy-power-sizing.json) uses
an explicit latent plateau, a two-stage heating pulse, a temperature limit,
ambient radiation, and the existing scalar workload search. It begins inside
the plateau and sizes the pulse against a limit above the melting temperature:

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/enthalpy-power-sizing.json
```

Its `transient.power_design` object declares `min_power_multiplier`,
`max_power_multiplier`, `power_multiplier_tolerance`, `temperature_tolerance_k`
and `max_evaluations`. A `sampled-peak` adjoint supplies the search with a checked
workload sensitivity. Remove `transient.adjoint` to use derivative-free sizing;
remove `transient.power_design` to evaluate just the declared schedule and its
gradient. Each candidate must run its complete coupled trajectory and satisfy
the physical acceptance gates before its sampled feasibility is used. The
returned feasible workload is checked across that trajectory's sampled
endpoints; it is not a continuous-time temperature certificate or a global
optimality claim. The example supplies synthetic numerical data, not a
validated hardware power limit.

The [enthalpy fan-sizing request](enthalpy-fan-sizing.json) instead declares
`hydraulics.fan`, an explicit `fan_speed_ratio` for each interval, and
`transient.fan_speed_design`. Its controls are `min_speed_multiplier`,
`max_speed_multiplier`, `speed_multiplier_tolerance`,
`temperature_tolerance_k` and `max_evaluations`; it cannot be combined with
`transient.power_design`.
The synthetic chart starts inside its latent plateau and the two-stage
workload increases from power scale 1 to 6, taking the endpoint into the liquid
sensible regime. Fan cooling during the early constant-temperature stages
changes stored enthalpy and therefore the later temperature: that history
contribution is retained in the fan gradient. The declared temperature limit
and speed bracket are illustrative numerical inputs, not binary-result or
hardware-performance claims.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/enthalpy-fan-sizing.json
```

Each candidate multiplies the complete declared fan schedule, starts from the
same initial enthalpy and recomputes actual airflow, convection and radiative
feedback over the whole trajectory. The forward solve, selected-speed report
and adjoint reconstruction all use those effective speeds. A sampled-peak
adjoint can suggest safeguarded search steps; missing smooth card derivatives
leave the speed slope unavailable and the search uses evaluated bisection
steps. Omitting `transient.adjoint` selects derivative-free sizing. The final
`transient_fan_speed_design` result reports the selected multiplier, effective
interval speeds, evaluated sampled peak and search history. It supplies no
electrical-power prediction, global minimum-speed proof or continuous peak
certificate. The fixture is synthetic and makes no measured fan-performance
claim.

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
