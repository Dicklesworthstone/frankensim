# Component allocation with independent temperature limits

```bash
frankensim --json cooling-component-design \
  examples/cooling-network/adjoint-component-contact-pulse.json \
  examples/cooling-network/allocate-multiple-temperature-limits.json
```

The allocator now requires EVERY declared temperature limit to pass on the
same selected power vector. The original base objective and allocation-level
`temperature_limit_k` remain the primary constraint. Optional extra constraints
can protect different components or surfaces at different temperatures.
A cooler memory component can be limiting even while a hotter chip passes its
own higher limit. The inputs and limits here are illustrative, not validated
hardware operating limits.

## Input and physical meaning

Add to the allocation document, not to the cooling base:

```json
"thermal_constraints": [
  {"name":"chip-limit","component":"chip","temperature_limit_k":302},
  {"name":"memory-limit","component":"memory","temperature_limit_k":301.4},
  {"name":"wall-limit","objective":{"mean_wall_region":"last-face"},
   "temperature_limit_k":301}
]
```

There can be 1..15 extra constraints, for at most 16 including the primary.
Names must be unique; `primary` is reserved. Each row requires exactly one
`component` or `objective`. A component selects the maximum over its original
PowerMap footprint vertices, including zero-power and overlapping components.
It does not infer an unmodeled junction temperature or expand the observation
region to every element touched by the heating source.

An explicit objective uses exactly one existing selector: `max_vertices`,
`max_solid_temperature: true`, `max_wall_region`, or `mean_wall_region`.
Unknown components/surfaces, empty or invalid vertex sets, repeated explicit
vertices, conflicting selectors, and gradient fields inside a selector refuse.
The existing cooling parser remains responsible for the complete physical input.
No constraints means the previous single-objective allocation path.

Each criterion observes its OWN peak over the original initial state and ALL
accepted endpoints of the fixed trajectory, including all repeated cycles.
Peaks can occur at different vertices and times. Neither the primary peak's
location/time nor the final temperature is substituted for another criterion.
Durations, step counts, initial temperatures, materials, heating footprints,
contacts, fan schedules and either radiation model remain unchanged.

## Search and derivatives

At a candidate, the search uses

```
maximum_temperature_excess_k = max_i(sampled_peak_i - temperature_limit_i).
```

Passing means this signed excess is nonpositive. It is not an average of
violations and not the largest absolute temperature. Both the original watt
bracket width and the slack of the most restrictive limit must meet their
respective stopping tolerances, unless the declared upper watt bound passes.
Priority ordering and conditioning on earlier selected/later minimum watts
retain the semantics described in COMPONENT_POWER_ALLOCATION.md.

Requested component adjoints are recomputed for EACH observed functional.
The active excess uses that criterion's absolute K/W slope, including its own
peak time and complete preceding thermal history. Exact ties between constraint
margins withhold a single slope and use bisection. Near ties and each underlying
spatial/time maximum retain the producer's selected-branch limitations.
Newton proposals remain inside the central 80% of the evaluated bracket; every
candidate is checked against every constraint. Predictions never decide passing.
No adjoint request means no hidden reverse solves.

## Cost, limits and interrupted work

This implementation runs the existing `cooling-network` producer once per
criterion per candidate. It is NOT a batched multi-objective forward solver.
The extra work is explicit: three criteria require three complete trajectories,
plus three criterion-specific adjoints when requested. The final physical fields,
transport quantities and integrated energy summaries must replay identically
across those observation changes; reverse-work counts and peak locations may differ.
No alternate thermal simulator or frozen-physics gradient is introduced.

`max_evaluations` counts candidate SETS, just as it previously counted candidate
vectors. `max_total_steps` counts every completed physical trajectory across all
criteria and priorities. The full next set is admitted against that cap before
launch. `wall_seconds` is one parent deadline for the complete allocation; every
child's original solver, timestep and adjoint-memory settings remain in force.
`trajectories_per_candidate` makes the extra cost visible.

`evaluations_attempted/completed` count candidate sets. The separate
`trajectory_evaluations_attempted/completed` count physical calls. Completed
trajectory steps and solid solves INCLUDE finished criteria of an interrupted
set. Work inside a killed, unfinished child is not invented. An incomplete set
never replaces a fully evaluated passing allocation. Parent-budget exhaustion
retains the previous passing vector with exit 6; a model failure, malformed
result or failed requested adjoint aborts without a design result.

Secondary full fields are dropped before the next criterion. Trial history
retains powers, primary peak, worst excess and active constraint; it does not
replicate large vertex selectors or all gradients on every trial. Full selected
observations are retained once. These bounds are not a total-process memory
certificate, and this feature adds no durable optimizer checkpoint.

## Results and replay

`cooling_result`, `selected_sampled_peak_k`, and `resolved_request` keep the
PRIMARY objective and its actual temperature. No shifted artificial Kelvin
value is published as a physical peak. `thermal_constraints` contains the
all-passed flag, maximum signed excess, active name, exact active count, and
one row per criterion. Rows retain their objective, limit, peak, peak time,
signed excess and derivatives in allocation-priority order.

To replay a secondary result, copy `resolved_request`, replace its `objective`
with that row's `objective`, and set `transient.temperature_limit_k` to that
row's limit. Run ordinary `cooling-network`. The temperature-limit field controls
reporting, not physics. The replay must match that criterion's peak/derivatives
and the primary result's final physical field. The allocation input files are
never rewritten.

This is deterministic nominal allocation on the declared fixed grid/cycle count,
not joint probabilistic reliability, a global/lexicographic optimum, a monotonicity
proof, mesh/time error control, physical validation, native `.fsim` integration,
or a continuous-time temperature certificate. Adaptive/controller/periodic-stop
and nested design/time/mesh-study restrictions of the original allocator remain.

## Focused verification

```bash
cargo test -p fs-cli --bin frankensim uq_command::component_design
cargo test -p fs-cli --test cooling_component_design
```

Fourteen new Rust tests comprise nine focused admission/search/replay tests and
five actual-command tests. They cover colder binding components, active-limit
switches and exact ties, independent peak times, initial wall maxima, each
criterion's own adjoint, unchanged physical fields, derivative-free execution,
closed-enclosure reflection, complete-set resource admission, interrupted-set
accounting and refusal of missing/bad constraints. Earlier single-objective
regressions remain. These Rust tests have NOT been executed here: the authoring
environment lacks a Rust toolchain. Compilation and actual CLI behavior remain
unverified. No new dependencies were added; the prior enclosure-adjoint lockfile
regeneration remains unresolved.

Independent NumPy/SciPy P1 calculations ran eight allocation cases covering
single/repeated trajectories, constant/nonlinear materials, matching/nonmatching
contact, and reservoir-radiating/nonradiating models. Separately solving a Brent
root for each constraint and intersecting those conditional intervals differed
from the chosen power by at most 0.000521 W, below the declared 0.001 W tolerance.
Thirty-six complete-trajectory derivative comparisons differed by at most
2.42e-11 K/W. These are independent mathematical references, NOT Rust execution.

The example selects standby=3 W, memory=6 W and chip=10.123629 W. The primary
and chip peaks are 301.781650 K; memory reaches 301.400000 K and is limiting.
Checking only the 304 K primary limit accepts chip=24 W, with a global peak of
303.861436 K but a 301.589985 K memory peak: both extra limits would fail.
Tightening the chip limit to 301.5 K instead selects about 8.259823 W and moves
the binding constraint to the chip. A separate heater-only calculation peaks
at the chip at 20 seconds, at memory at 28 seconds, and at the observed mean
wall initially. The guided example used 18 candidate sets, corresponding to
54 producer trajectories and 756 accepted endpoints in this implementation;
these are reference operation counts, not measured Rust performance.
