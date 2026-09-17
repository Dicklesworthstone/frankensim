# Full-trajectory timestep studies

```bash
frankensim --json cooling-network \
  examples/cooling-network/time-convergence-radiative-pulse.json
```

`transient.time_convergence` runs the actual cooling trajectory repeatedly on
nested fixed time grids. Each level doubles the number of steps in EVERY
workload/fan interval, retaining interval durations and switch times. It compares
all nodal temperatures at every common endpoint, AND the selected objective's
sampled peak over the entire trajectory, including newly introduced endpoints.
Both changes must be within the declared Kelvin tolerance for at least two
consecutive comparisons. A cooled final temperature or an unchanged initial
peak cannot stand in for the complete temperature history.

Add inside an ordinary `transient` request:

```json
"time_convergence": {
  "max_refinements": 6,
  "consecutive_passes": 2,
  "temperature_tolerance_k": 0.05,
  "max_total_steps": 5000,
  "max_trace_bytes": 16777216
}
```

The example uses two 150-second pulse/cooldown cycles, nonlinear material
conductivities, heterogeneous heat capacities, finite contact resistance,
radiation, and fan-driven bypass mixing. Its data are illustrative declarations,
not measured hardware. The 0.05 K setting is an observed-comparison policy, not
an error certificate. A run can legitimately exhaust its original wall allowance.

## The same physical experiment on every grid

Every level starts with the ORIGINAL nodal initial temperatures. Only cycles
within that level inherit heat from earlier cycles. Component footprints, named
interval powers, material laws, contact operators, emissivities, surroundings,
airflow controls and initial data stay unchanged. The existing endpoint producer
solves storage and all declared nonlinear boundary physics. Its accepted-step,
whole-window and repeated-cycle energy checks still apply.

A read-only observer receives accepted temperatures after their physical checks.
It cannot change the previous field or record a Newton, radiation or air-coupling
trial as a timestep. Global observation time is used across repeated cycles;
physical integration retains the existing local cycle time. No interpolation
crosses a workload or fan switch. Common timestamps must agree bit for bit.

The final published `transient`/`repeated_cycles` result is the finest completed
passing trajectory only. Its input, stored, radiative and exhaust energy are NOT
summed across the exploratory grids. `time_convergence.total_accepted_steps` and
`total_solid_solves` separately count all completed trajectory work.

Matching and explicitly declared nonmatching planar contacts use their existing
operators. This feature does not refine either volume mesh or contact geometry.
It is a temporal study, not simultaneous space-time adaptation.

## Limits are cumulative and explicit

`transient.max_steps` still caps accepted steps per cycle. A declared
`repeat.max_total_steps` still caps one complete repeated trajectory. The NEW
`time_convergence.max_total_steps` caps their sum across all solved grids; it
does not replace or increase either original cap. The next grid's complete step
count is checked before it runs. `max_refinements` limits subdivisions after the
base trajectory. Too few successful comparisons produce a budget refusal with
no convergence result, not permission to loosen the tolerance.

`max_trace_bytes` caps the two retained accepted-field buffers, including their
row times and buffer headers. It is NOT a total-process memory limit: the normal
solver, schedule and result buffers remain separately owned. Allocation uses
checked arithmetic and fallible reservation. Only adjacent grid traces remain
in memory. All trajectory solves share the original command wall deadline.
No new per-grid time allowance or inner iteration budget is granted.

This mode admits only fixed schedules and fixed cycle counts. Adaptive grids,
periodic stopping, hysteretic controllers, adjoints and nested nominal sizing
with this policy refuse explicitly. Their ordinary non-study paths are not
changed. A future derivative of an automatically chosen grid would need its
own semantics; this feature does not pretend to provide one.

## Replay the reported finest grid

The result supplies `time_convergence.final_steps_per_interval`. Remove
`transient.time_convergence` and assign those integers to each corresponding
`transient.intervals[i].steps`, leaving all other input unchanged. The normal
command then replays the final trajectory without the ladder or trace buffers.

`interval.steps` is optional. Without it, the previous rule
`ceil(duration_s / max_step_s)` is unchanged. An explicit count must be an integer,
positive, and at least that minimum. For example, duration 2.3 seconds with
max_step_s 1 starts with three steps. Its refinement has SIX, not the five that
would result from simply halving max_step_s and repeating the ceiling operation.
This distinction preserves common endpoints and supports exact final-grid replay.

## Meaning of success

The status is `successive-time-grid-tolerance-met`. History records interval
counts, all-cycle peak and time, common-endpoint field differences, peak changes,
work and consecutive passes. Newly introduced endpoints affect the peak but
have no coarse full-field partner; the next refinement subsequently checks them.
Peak times may move. No monotonic convergence, Richardson error bound,
continuous-time maximum or threshold-crossing certificate, spatial accuracy,
physical validation or hardware compliance is inferred from agreement.

## Focused checks

```bash
cargo test -p fs-cli --test cooling_time_convergence
cargo test -p fs-cli --bin frankensim time_convergence
cargo test -p fs-cli --test cooling_transient_radiation
```

Ten new Rust regressions comprise four direct policy/grid tests and six actual
command tests. They cover complete field comparison when maxima agree, exact
common-time alignment, noninteger interval lengths, repeated radiating storage,
nonmatching contact, final-grid replay, default-path equivalence, cumulative and
per-trajectory limits, and unsupported policies. These Rust tests were authored
but NOT executed in the authoring environment, which has no Rust/RCH toolchain.

Independent NumPy P1 calculations with direct nonlinear solves and analytic
air elimination exercised the same illustrative radiating model. For the full
two-cycle example, the base 30-endpoint trajectory peaked at 306.792379 K. At
level 3 its peak change was 0.037066 K but its full-field change was 0.066568 K:
the peak alone would pass 0.05 K too early. Levels 4 and 5 passed both tests;
the final 960-endpoint trajectory peaked at 307.070667 K at 180 seconds, with
0.020148 K common-field change and 0.009422 K peak change. All six trajectories
performed 1,890 accepted endpoints in this independent calculation.

A separate hot-initial mean-wall case had unchanged 350 K sampled peaks but
0.077093 K full-field change. Ten constant-material controls compared backward
Euler against independent matrix-exponential solutions of the SAME fixed-mesh
semidiscrete system; all tested refinements reduced maximum temporal field
error. For the contact control it decreased from 0.354715 to 0.040170 K over
five grids. These are independent numerical references, not executions of
FrankenSim, continuum-space bounds or experimentally validated temperatures.
