# Size the fan schedule against a workload trajectory

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/size-transient-fan.json
cargo test -p fs-cli --bin frankensim network_command::transient::sizing
```

The existing steady `fan_speed_design` cannot answer whether a changing workload
stays below a temperature limit before it cools down. A **nested**
`transient.fan_speed_design` searches against the peak sampled objective over
the complete transient. Its objective can be a surface mean, surface maximum,
whole-solid maximum, or selected-node maximum, as in the existing command.

Add this object inside the existing `transient` request:

```json
"temperature_limit_k": 320.85,
"fan_speed_design": {
  "min_speed_multiplier": 0.5,
  "max_speed_multiplier": 1.3,
  "speed_multiplier_tolerance": 0.0001,
  "temperature_tolerance_k": 0.00001,
  "max_evaluations": 64
}
```

The one selected multiplier scales **every interval's declared fan speed**. It
does not replace the schedule with a constant speed. Interval durations,
component powers, thermal capacities, geometry, and initial temperatures stay
unchanged. Both bound multipliers must keep every interval inside the fan's
speed domain; this is checked before physics. Each actual trial still has to
pass fan-curve, hydraulic and convection validity checks. No regime switching,
coefficient clipping or extrapolation is introduced.

Every candidate starts from the same initial field and runs the actual fan,
flow-derived convection, contact-aware solid and mixed-air trajectory. No state
from a previous candidate leaks into the next. Feasibility uses the initial
objective and all accepted endpoint samples, including both half-step samples
under adaptive integration. Neither the final field nor a steady-state
surrogate decides whether a trial passes.

The bracket begins with an evaluated failing lower multiplier and passing upper
multiplier. By default, bisection retains the passing trial, returning it only
when both the multiplier width and passing-side temperature slack meet their
tolerances. A passing declared minimum returns after one complete evaluation.
If neither endpoint passes, the result is a missing bracket, not a proof that
every interior speed fails. The result is not a global minimum-speed certificate.

The output contains the passing trial's complete transient history and final
field, including matching final fan flow and coefficients. Its
`transient_fan_speed_design` object adds the multiplier, actual interval speeds,
failed lower endpoint, evaluated peak/time history, and work totals across all
candidates. Only the passing trajectory is retained; other trial histories are
scalar summaries. Steady adjoints remain disabled. Transient peak adjoints are
separately opt-in as described below; no electrical-energy claim is implied.

`max_evaluations` is a whole-trajectory budget (at most 256). Existing step and
adaptive-trial limits apply separately to each candidate. The original single
wall budget covers the entire design search. A numerical/domain/cancellation
failure aborts instead of being treated as a failing temperature candidate.
Evaluation exhaustion publishes no partial trajectory and uses exit class 6.

## Optional sampled-peak adjoint guidance

A fixed-grid, fixed-cycle-count design may add an explicit
`"adjoint":{"qoi":"sampled-peak","max_checkpoint_bytes":1048576}` beside its
fan or workload design object. Keep `objective.gradient` false. This enables
candidate-specific discrete trajectory adjoints and guarded Newton proposals,
with bisection fallback. It does not change either acceptance tolerance, allow
predicted feasibility, or silently increase the resource budgets.

```bash
frankensim --json cooling-network examples/cooling-network/size-transient-adjoint-fan.json
frankensim --json cooling-network examples/cooling-network/size-transient-adjoint-power.json
```

Every candidate binds its actual schedule to both the forward simulation and
reverse reconstruction. Repeated warm-up history contributes to the gradient
and to feasibility. `search_method`, `newton_trials`, and each history row's
`dpeak_dmultiplier_k` expose the policy and local outer-coordinate slope.
Without the explicit adjoint policy, bisection remains and no reverse work is
performed. Guidance is not guaranteed to reduce evaluation count or runtime.

Final-state adjoints, adaptive grids, periodic stopping and hysteretic controller
adjoints cannot be used for this peak design path. Ordinary derivative-free
runs retain those existing simulation capabilities. Requested adjoint failures
abort; an unavailable or unsuitable proposal slope falls back to bisection.
See [TRANSIENT_ADJOINT.md](TRANSIENT_ADJOINT.md) for candidate-relative versus
outer-multiplier units, repeated controls, memory accounting and exact replay.

## Time resolution changes the engineering decision

The `size-transient-fan.json` sample is a two-material, contact-coupled slab
heated locally with 20 W for 30 seconds, then cooled for 120 seconds. It uses
the existing duct correlations and a 1 / 1.5 base fan schedule. Illustrative
heat capacities are 200,000 and 100,000 J/(m3 K); they are declarations, not
material-card claims.

Independent NumPy P1/contact/air calculations give approximately:

| Integration | Baseline sampled peak | Selected fan multiplier |
|---|---:|---:|
| Fixed 2-second steps | 320.736502 K | 0.686035 |
| Adaptive step doubling as declared | 320.884207 K | 1.108398 |

For the 320.85 K target, coarse fixed steps accept the baseline while the
adaptive run rejects it. The committed request uses the adaptive settings.
Its reference passing peak is 320.849991 K at 30 seconds; the final maximum is
only 302.887171 K. Checking just the cooled final state would miss the problem.
These are independent mathematical references, **not Rust execution results**.

Adaptive step doubling estimates local endpoint error. It does not certify
midpoint errors, the complete trajectory, or extrema between accepted samples.
Different candidates may use different adaptive sample times, so a narrow fan
bracket is not a physical compliance margin or a timestep-convergence proof.
Results remain nominal discrete estimates, with quasi-steady air, frozen fluid
properties, fixed geometry/contact resistance, no fan-off natural convection,
and no native `.fsim`/ledger integration. Rust compilation and tests were not
executed in the authoring environment.

## Size workload power while keeping the cooling schedule fixed

The complementary `transient.power_design` searches for a passing workload
multiplier near the temperature crossing, with an evaluated **failing upper**
endpoint. It is mutually exclusive with `transient.fan_speed_design` and with
both existing steady design modes.

```json
"temperature_limit_k": 315,
"power_design": {
  "min_power_multiplier": 0,
  "max_power_multiplier": 2,
  "power_multiplier_tolerance": 0.0001,
  "temperature_tolerance_k": 0.00001,
  "max_evaluations": 64
}
```

The common factor scales every interval's global `power_scale`, or every
absolute wattage in that interval's `component_powers_w` map. This preserves
the relative workload distribution, component footprints and durations. Zero
powers remain zero; there is no division by nominal component power. Uniform
signed source fields are scaled as declared, including cooling sources. Initial
temperatures, inlet temperatures, thermal capacities, and fan speeds do not
scale. Named watts are reprojected through the existing `PowerMap`, not guessed
from the already summed source. Pressure-driven schedules are supported too.

The search evaluates the requested maximum first. If it passes, it returns
`maximum-feasible`. Otherwise a passing lower endpoint is required. The search
returns the **passing lower** field and workload, never the failing upper trial.
Bisection is the default; the explicit sampled-peak adjoint policy can guide
interior trials under the same acceptance rules. Missing brackets and producer
failures retain their previous meanings; no global monotonicity, throughput
model, maximum safe hardware rating or continuum compliance certificate is
inferred. The same sampled-trajectory limitations and whole-search wall and
evaluation budgets apply.

`transient_power_design` reports the selected power multiplier, failed upper
endpoint, actual interval workload maps, peak/time and work totals. The main
transient history reports the applied watts or global scales, while
`solid_inputs` remains the base declaration. The final fan result corresponds
to the original unscaled final interval speed. Only the chosen passing
trajectory is retained in full.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/size-transient-power.json
```

This request uses explicit component watt maps, the same contact/duct/solid
model, and the fixed 1 / 1.5 fan schedule. For a 315 K sampled-peak target,
an independent adaptive NumPy reference selects a multiplier near 0.71838045:
the 20 W pulse becomes approximately 14.367609 W, with a passing peak near
314.999998 K and about 431.028271 J of input. Its cooled final maximum is
302.149715 K. These are independent mathematical checks, not Rust execution.
The added Rust tests for workload scaling and design remain unexecuted.
