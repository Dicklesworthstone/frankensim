# Repeated workload cycles without a cold reset

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/repeated-contact-pulse.json
cargo test -p fs-cli --bin frankensim network_command::transient::repeat
```

The optional `transient.repeat` object repeats the entire existing interval list:

```json
"repeat": {"cycles":10,"max_total_steps":10000}
```

`cycles` is an integer from 1 through 4096. The thermal field at the end of one
cycle is the initial field of the next. Capacities, contact traces, named source
footprints, interval lengths and fan schedules remain unchanged. Each cycle runs
the existing backward-Euler solid and quasi-steady air solve. No state is inferred
from a scalar hotspot and no heat is erased at a cycle boundary.

Both fixed and adaptive stepping are supported. Fixed timesteps preserve the
original per-interval partition. Adaptive stepping restarts its step suggestion
at each interval boundary, as before; only accepted fields and heat become
history. `max_steps` and adaptive `max_trials` remain per-cycle budgets.
`max_total_steps` adds a whole-run accepted-endpoint cap (at most 1,000,000).
That remaining budget is passed into the next cycle before it advances, not
checked after an unlimited cycle. Cycles and the original wall budget bound the
whole invocation. A refused solve or exhausted budget publishes no partial run.

The top-level final temperature field, fan result and coefficients describe the
last accepted cycle endpoint. To avoid storing every nodal trajectory, the
existing `transient` record contains only the last cycle's detailed history in
**local cycle time**. `repeated_cycles.last_cycle_start_time_s` supplies its global
offset. `repeated_cycles.cycles` contains scalar summaries of every cycle,
including its peak and global peak time, full-field start/end difference, heat
storage, input and exhaust. Whole-run totals count all cycles. Internal contact
transfer is not counted again as external heat.

`repeated_cycles.sampled_peak_objective_k` includes warm-up, the declared initial
condition and every accepted endpoint, not just the last cycle. The first
sampled limit violation likewise uses global elapsed time. A hot initial
condition or earlier overshoot cannot disappear from the reported maximum.
No limit on unsampled times, future cycles or the continuum field is inferred.

## Repeated-cycle sizing

The existing nested `transient.fan_speed_design` and `transient.power_design`
operate on the entire repeated experiment. Every design candidate starts from
the original declared initial field. Cycles inside that candidate inherit heat;
candidates do not inherit heat from each other. The power-scaling path retains
the repeat configuration, including with independent named component workloads.
Design summaries and work counts cover every cycle. A returned design is an
evaluated passing sampled trajectory, not a globally optimal or physically
certified operating limit. Parameter bounds and physical domain refusals are
unchanged. The total step cap applies per design candidate and the wall budget
covers the complete search.

## Example and verification status

The example repeats the existing 20 W, 30-second contact-coupled pulse followed
by 120 seconds of cooling, with a 1 / 1.5 fan schedule. An independent NumPy P1
FEM calculation on the declared mesh predicts a first-cycle peak of about
306.343159 K, but 310.748442 K in the tenth cycle, at 1380 seconds. A 308 K limit
passes the first pulse and fails in the third. The total input is 6000 J;
the final solid still reaches only about 304.948465 K after cooling.
FrankenSim itself (measured 2026-09-25, release build, after 05db922bf made the transient capacitance row-sum lumped; the reference above uses the consistent P1 mass, so the two differ on this 12-tet mesh) gives a first-cycle peak of 304.1166 K
and 309.1035042 K in the tenth cycle, at 1380 seconds. The example's limit is
now 306 K: the first two cycles pass (304.12 and 305.63 K), and the third first
exceeds it at 326 seconds. The final solid reaches 305.5585429 K.

These are independent mathematical references, not executed Rust measurements.
The focused Rust tests cover state carryover, energy totals, single-cycle parity,
power sizing and refusal behavior. Rust compilation, formatting and test execution
were unavailable in the authoring environment. This remains the experimental
JSON cooling command, not a native `.fsim`/ledger workflow. Repeating a specified
number of cycles does not establish a settled periodic thermal state.

## Run until successive cycles agree at the same phase

Instead of `cycles`, supply `until_periodic`:

```json
"repeat": {
  "until_periodic": {
    "max_cycles":100,
    "temperature_tolerance_k":0.0001,
    "consecutive_cycles":2
  },
  "max_total_steps":10000
}
```

The condition is `max_i |T_end[i] - T_start[i]| <= temperature_tolerance_k`
over **every solid node**, checked at the same phase of each complete schedule.
An unchanged surface mean or maximum does not pass this test when heat is still
redistributing internally. The condition must pass on at least two consecutive
cycles; a failure resets the streak. No relaxation or extrapolation alters the
accepted field, and no endpoint from a partially completed cycle can pass.

`max_cycles` is a hard cap from 2 through 4096. `consecutive_cycles` must be at
least 2 and no greater than that cap. The total step cap must accommodate the
minimum number of qualifying cycles; it can be smaller than the worst-case
`max_cycles * max_steps`. Exhausting either budget returns a refusal rather than
publishing an unconverged cycle as periodic. Solver and cancellation failures
propagate through the same boundary. Capacities, source projection, contact,
fan curves and correlation domains use the original producers on every cycle.

The result status is `periodic-field-tolerance-met`. The `periodic` object retains
the actual maximum nodal residual, tolerance and achieved consecutive count.
`cycles_completed` is the number actually executed, not the declared cap.
The small final-cycle stored energy and corresponding heat imbalance remain
reported as computed: they are not forced to zero to manufacture periodicity.
The full warm-up peak remains the design objective even when the final cycle
is cooler, and every design candidate independently meets the periodic gate.

**This is an observed residual of the discrete cycle map, not an error bound on
the infinite-cycle solution.** A slowly contracting mode can have a small change
per cycle while remaining farther from its limit. Adaptive sample times can
also change between cycles. The gate is not a waveform-error bound, proof of
uniqueness or stability, guarantee on later-cycle maxima, or physical temperature
margin. Tightening it does not replace timestep or spatial refinement. A
periodicity failure is not evidence that a proposed power or fan setting is
thermally infeasible; sizing propagates that failure.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/periodic-contact-pulse.json
cargo run -p fs-cli --bin frankensim -- --json cooling-network \
  examples/cooling-network/size-periodic-power.json
```

The independent fixed-step P1 reference meets the 0.0001 K full-field criterion
on cycles 51 and 52. The accepted cycle peaks near 311.672519 K. Direct solution
of the independently assembled affine cycle fixed point gives approximately
311.672953 K: even in this benign example, the distance to the limiting waveform
is not the same as the accepted cycle-map residual. Power sizing against a
308 K limit selects a multiplier near 0.6853, reducing the 20 W pulse to about
13.71 W. A first cold pulse alone would have accepted the original 20 W.
These numbers remain independent mathematical references, not Rust executions.
FrankenSim with lumped capacity (measured 2026-09-25) accepts cycle peaks near
310.0425874 K, and power sizing against 308 K selects a multiplier near 0.79657,
about 15.93 W.
