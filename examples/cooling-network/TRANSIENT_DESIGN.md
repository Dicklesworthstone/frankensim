# Size the fan schedule against a workload trajectory

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/size-transient-fan.json
cargo test -p fs-cli --bin frankensim network_command::transient::sizing
```

The existing steady `fan_speed_design` cannot answer whether a changing workload
stays below a temperature limit before it cools down. A new **nested**
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
multiplier. Bisection retains the passing trial, returning it only when both
the multiplier width and passing-side temperature slack meet their tolerances.
A passing declared minimum returns after one complete evaluation. If neither
endpoint passes, the result is a missing bracket, not a proof that every
interior speed fails. The result is not a global minimum-speed certificate.

The output contains the passing trial's complete transient history and final
field, including matching final fan flow and coefficients. Its
`transient_fan_speed_design` object adds the multiplier, actual interval speeds,
failed lower endpoint, evaluated peak/time history, and work totals across all
candidates. Only the passing trajectory is retained; other trial histories are
scalar summaries. Steady adjoints remain disabled, and no fan-speed gradient or
electrical-energy claim is implied.

`max_evaluations` is a whole-trajectory budget (at most 256). Existing step and
adaptive-trial limits apply separately to each candidate. The original single
wall budget covers the entire design search. A numerical/domain/cancellation
failure aborts instead of being treated as a failing temperature candidate.
Evaluation exhaustion publishes no partial trajectory and uses exit class 6.

## Time resolution changes the engineering decision

The sample is a two-material, contact-coupled slab heated locally with 20 W
for 30 seconds, then cooled for 120 seconds. It uses the existing duct
correlations and a 1 / 1.5 base fan schedule. Illustrative heat capacities are
200,000 and 100,000 J/(m3 K); they are declarations, not material-card claims.

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
