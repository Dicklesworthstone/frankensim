# Cycle-rate fan feedback for repeated thermal workloads

Repeated transient cooling can optionally use a hysteretic fan controller driven
by one explicitly declared solid-mesh sensor vertex. This is a closed-loop
extension of `transient.repeat`: the accepted end field of one cycle is both the
thermal history and the sensor state used to select fan speed for the next
cycle.

Add `fan_controller` inside `transient.repeat`:

```json
"repeat": {
  "cycles": 10,
  "max_total_steps": 10000,
  "fan_controller": {
    "sensor_vertex": 4,
    "low_temperature_k": 300.2,
    "high_temperature_k": 301.0,
    "low_speed_multiplier": 0.7,
    "high_speed_multiplier": 1.3,
    "initial_speed_multiplier": 0.7
  }
}
```

The controller requires `hydraulics.fan`. The selected multiplier scales every
fan speed already declared by the cycle's interval schedule; it does not erase
relative interval speeds. All low, high and initial combinations are checked
against the fan's declared speed domain before the first controlled cycle, and
the actual fan/graph and flow-derived convection producers still run normally.

At the start of each cycle:

* sensor >= `high_temperature_k` selects the high multiplier;
* sensor <= `low_temperature_k` selects the low multiplier;
* between thresholds, the previous multiplier is retained.

The selected multiplier is held for the complete cycle. There is no sub-cycle
sampling, sensor filtering, noise, delay, actuator inertia, PWM model, inferred
RPM, electrical-power estimate, or hidden interpolation. The first cycle begins
from `initial_speed_multiplier`, but its declared initial temperature is sampled
before work and may switch it immediately.

Each cycle summary reports start/end sensor temperatures, controller and total
fan multipliers, whether a switch occurred, and whether the end state would keep
the same hysteretic controller state on the following cycle. The aggregate
`fan_controller` result reports the threshold/speed declarations and switch
count.

When `repeat.until_periodic` is used, convergence requires BOTH the existing
full-nodal same-phase temperature residual and an unchanged controller state.
A field that meets the temperature tolerance but would cross a thermostat
threshold on the next cycle is not declared periodic.

The outer `transient.fan_speed_design`, when present, remains an overall scale
on top of the controller-selected multiplier; every design candidate starts
from the original initial field and fresh controller state. Domain failures are
refusals, not failing temperature candidates. `transient.power_design` can
likewise evaluate a fixed controller while scaling the workload.

This controller operates only at complete duty-cycle boundaries. It is useful
for cycle-rate supervisory cooling policies and periodic workload studies, not
for millisecond-scale thermal control. The existing sampled-temperature,
quasi-steady-air, timestep, continuum, and physical-validation limitations
remain unchanged.
