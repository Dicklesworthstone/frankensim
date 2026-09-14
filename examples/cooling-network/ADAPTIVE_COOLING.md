# Adaptive timesteps for coupled cooling

The existing `cooling-network` transient can now choose timesteps using the
**full solid-temperature field**, not only the selected hotspot or surface mean.
It retains the same heterogeneous FEM, finite contact resistance, component
workloads, fan schedule and quasi-steady mixed-air transport.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/adaptive-contact-pulse.json
cargo test -p fs-cli --bin frankensim network_command::transient::adaptive
```

Inside the existing `transient` object, add:

```json
"adaptive": {
  "absolute_tolerance_k": 0.01,
  "relative_tolerance": 0.001,
  "minimum_trial_step_s": 0.00001,
  "max_trials": 4000
}
```

All four fields are required. The absolute tolerance and minimum trial width
must be positive. The relative tolerance is in `[0,1)` and the trial count is
in `1..=100000`. The minimum trial width cannot exceed `max_step_s`. Omission of
`adaptive` preserves the existing fixed-step trajectory and numerical ordering;
its result's `adaptive` member is `null`.

## What is evaluated and accepted

For an attempted interval of width H, the controller evaluates one full
backward-Euler step and two successive half-steps from the same accepted
previous field. Each evaluation fully converges the solid/air coupling using
that evaluation's immutable previous field. Source, geometry, properties,
contact law, fan flow and convection are constant across the trial; no trial
crosses a workload or speed discontinuity.

For each vertex i, define

```
change_i = max(abs(T_coarse_i - T_old_i), abs(T_fine_i - T_old_i))
ratio_i = abs(T_fine_i - T_coarse_i) / (absolute_tolerance_k + relative_tolerance * change_i)
error_ratio = max_i ratio_i
```

Relative scaling uses temperature **changes**, not a large absolute-kelvin
offset. The first-order step-doubling discrepancy estimates the fine endpoint's
local error. When `error_ratio <= 1`, the two half-step results are accepted
without Richardson extrapolation. Both the midpoint and final endpoint enter
sampled hotspot/limit detection. The error estimate belongs to the **final**
endpoint, not to an independent test of the midpoint.

Otherwise, both trial fields and their heat integrals are discarded, H shrinks,
and the same old field is retried. The controller uses a safety factor of 0.9
and square-root error scaling, limits growth to 2, and limits rejected-trial
shrink factors to [0.2,0.8]. It resets its width proposal at each input interval.
A solver, energy, correlation-domain, cancellation or material failure is not
reinterpreted as a truncation-error rejection: that producer failure propagates.

`max_step_s` and `minimum_trial_step_s` govern the **full trial** width; each
accepted integration substep is half that width, up to floating-point endpoint
rounding. A final interval-clipped trial may be shorter than the minimum. If it
still fails the error test, the run refuses rather than stepping below the
minimum to conceal failure. Unrepresentable midpoints also refuse.

## Budgets and accounting

`max_trials` counts every full-versus-two-half comparison, including rejected
trials, across the whole invocation. `max_steps` counts accepted half-step
endpoints; the complete numerical wall budget still covers all work. A trial
requires room for two endpoints before running. Exhaustion publishes no partial
trajectory. The initial state does not consume an integration step.

`total_solid_solves` includes discarded coarse and rejected solves. Heat storage,
input energy and external-air heat include **only the accepted half-steps**.
Contact heat remains internal and is not counted again. Result history contains
`estimated_local_error_ratio` at each accepted pair's final endpoint and `null`
at its midpoint. The `adaptive` summary contains attempted, rejected and accepted
trial counts, tolerances, and the largest accepted error ratio.

## Scope and reference

This is **local error estimation**, not a rigorous temporal bound, a cumulative
global-error guarantee, a midpoint-error bound, or a continuum-mesh certificate.
Sampled peaks and first recorded violations still do not bound between-step
behavior or locate exact crossing times. Tightening the tolerance and checking
mesh/time refinement remains a resolution study, not physical validation.

The committed example is the existing 20 W contact-coupled pulse followed by
cooldown, with a 30-second maximum trial width instead of manual two-second
steps. Independent NumPy FEM calculations at the declared adaptive tolerances
produce 76 accepted half-step endpoints, 46 trials (8 rejected), a sampled peak
of approximately 306.368112 K at 30 seconds, and a final maximum of 301.485965 K.
The 600 J input splits into about 536.448925 J stored and 63.551075 J exhausted.
These are independent mathematical references, **not executed Rust results**.

Compilation, Rust tests and formatting have not been run in the authoring
environment. No new dependency, native `.fsim` integration, ledger workflow,
transient adjoint or fan-off convection model is introduced.
