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

These are independent mathematical references, not executed Rust measurements.
The focused Rust tests cover state carryover, energy totals, single-cycle parity,
power sizing and refusal behavior. Rust compilation, formatting and test execution
were unavailable in the authoring environment. This remains the experimental
JSON cooling command, not a native `.fsim`/ledger workflow. Repeating a specified
number of cycles does not establish a settled periodic thermal state.
