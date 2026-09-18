# Priority-ordered component power allocation

```bash
frankensim --json cooling-component-design \
  examples/cooling-network/adjoint-component-contact-pulse.json \
  examples/cooling-network/allocate-component-pulse.json
```

This command allocates independent component watts against the selected
objective's peak over a COMPLETE cooling trajectory. It does not multiply all
components together. The example activates a dormant standby footprint, gives
memory its declared priority, and uses the remaining sampled-temperature margin
for chip power. It reuses the existing nonlinear, contact, radiation and mixed-air
cooling producer. Inputs and allocation limits are illustrative declarations,
not measured hardware or validated operating limits.

## Declare the allocation policy

```json
{
  "schema": "frankensim.cooling-component-design.v1",
  "units": "SI",
  "temperature_limit_k": 302.5,
  "power_tolerance_w": 0.001,
  "temperature_tolerance_k": 0.0005,
  "max_evaluations": 128,
  "max_total_steps": 5000,
  "wall_seconds": 900,
  "priority": [
    {"component":"standby","interval":0,"min_power_w":0,"max_power_w":3},
    {"component":"memory","interval":0,"min_power_w":0,"max_power_w":10},
    {"component":"chip","interval":0,"min_power_w":0,"max_power_w":24}
  ]
}
```

Each control names one original `solid.component_power` footprint and one
zero-based schedule interval. Its absolute watts apply to EVERY occurrence of
that interval in a fixed repeated schedule. The same component can be controlled
in different intervals, but repeating the same component/interval pair refuses.
There must be 1..64 controls with finite bounds `0 <= min_power_w < max_power_w`.
The original component can have zero watts; no ratio to its original power is
needed, and overlapping footprints remain independently controllable.

Only controlled intervals are expanded to `component_powers_w`. Every unlisted
component retains its actual watts: a `power_scale` interval uses scale times
base watts; an absolute workload retains its explicit overrides. Other intervals
remain unchanged. Geometry, source footprints, materials, contacts, radiation,
fan settings, durations, initial temperatures and solver settings are preserved.
The resolved request's `transient.temperature_limit_k` is set to this allocation
limit for consistent threshold reporting; that field does not change the PDE.

The first candidate puts all controlled axes at their declared minima. If it
fails the temperature limit, the command refuses; this is NOT a proof that no
other allocation is feasible. Otherwise axes are visited in declaration order.
Earlier axes keep their selected values, and later axes remain at their minima.
An evaluated passing upper bound is accepted directly. Otherwise the current
passing allocation and failed upper trial define a conditional scalar bracket.
Both its watt width and the passing temperature slack must meet their tolerances.

This is a local, declared-order allocation policy, not a general constrained
optimizer or a proof of lexicographic/global optimality. Later allocations can
change the earlier conditional brackets. Every final selected VECTOR is still
an actually evaluated passing trajectory, never a stitched combination of
independently passing scalar results. The receipt explicitly describes each
priority decision's conditioning. Changing priority order can change the answer.

## Optional total-adjoint guidance

The base must have `objective.gradient=false`. With no `transient.adjoint`, the
search uses bisection and performs no hidden reverse solves. To guide candidate
selection, request `transient.adjoint` with `qoi:"sampled-peak"`,
`component_power:true`, and the existing explicit `max_checkpoint_bytes` budget.
The example base already contains this declaration. Final-temperature adjoints
cannot guide an all-cycle peak constraint and are rejected.

The driver verifies that the returned adjoint describes the SAME peak value,
peak time, cycle count, interval, component and applied watts as the candidate.
Its absolute K/W slope remains useful at zero watts. Finite positive slopes can
propose Newton trials only strictly within the central 80% of the evaluated
bracket; unusable slopes fall back to bisection. The derivative of a selected
maximum branch does not certify branch stability or finite changes. Every trial
runs the full physical trajectory, including the requested adjoint. Missing or
failed adjoints are errors, not permission to publish a frozen-physics slope.

## Resource limits and useful partial results

`max_evaluations` (1..4096), `max_total_steps` (1..100000000), and `wall_seconds`
(positive, at most 86400) apply across ALL priorities, not separately to each.
The next complete fixed trajectory is admitted against the cumulative step cap
before launch. Each child's original per-cycle, repeated-cycle, solver and
checkpoint-memory limits also remain in force. The parent deadline starts after
file parsing/allocation admission and covers candidate serialization, child
execution and result preparation. The existing same-executable child watchdog
kills and reaps interrupted children and joins their pipe workers.

Completed searches return `priority-allocation-complete` with exit 0. Allocation
budget exhaustion after a passing baseline returns `budget-exhausted` with exit 6
AND the last complete passing allocation. Watt representability exhaustion also
returns this noncomplete status with its reason. Before a passing baseline exists,
there is no partial result to publish. Model failures, malformed child results,
and a child's own numerical/solver-budget refusal abort without a design result;
they are never classified as hot samples or silently skipped.

`evaluations_attempted` includes an interrupted child; `evaluations_completed`,
`total_completed_trajectory_steps`, and `total_completed_trajectory_solid_solves`
count only complete trajectories. They do not invent a work count for a killed
child. Output formatting cannot upgrade an expired allocation deadline to success.
The input files are never rewritten. Controlled-interval expansion is capped at
65536 component rows; this is not a total-process memory guarantee.

## Replay and supported scope

Both complete and useful partial outputs retain `resolved_request` and
`cooling_result`. Save the former as a JSON file and pass it to ordinary
`cooling-network` to replay the exact selected trajectory, without rerunning the
allocation search. The original base source declarations remain in that request;
its resolved interval workloads are the actual selected watts. This is NOT a
durable optimizer checkpoint or mid-trajectory resume facility.

Fixed timesteps, explicit interval step counts, single cycles, fixed repetition,
nonlinear materials, matching/nonmatching contacts, and reservoir/enclosure
radiation use the same underlying producer. Nested searches, adaptive timesteps,
automatic time/mesh studies, periodic stopping and thermostatic controllers
refuse here rather than acquiring undocumented optimization semantics.
The command currently requires Unix for the existing cooling child interface.
It does not add probabilistic reliability, continuous-time peak certification,
mesh/time error bounds, physical validation or native `.fsim` ledger integration.
Sharing the UQ child runner does not make this a statistical allocation method.

## Focused verification

```bash
cargo test -p fs-cli --bin frankensim uq_command::component_design
cargo test -p fs-cli --bin frankensim uq_command::child
cargo test -p fs-cli --test cooling_component_design
```

Fourteen new Rust regressions comprise eight allocation-unit tests, one shared
child-output test, and five actual-command tests. They cover workload preservation,
zero-power controls, priority order, global budgets, wrong/missing adjoints,
selected-request equality, complete replay, enclosure reflection and useful
partial results. These tests have NOT run in the authoring environment, which
has no Rust toolchain. Compilation and actual CLI behavior remain unverified.
This feature adds no new dependencies; it does not resolve the previously noted
unregenerated lockfile after the enclosure-adjoint dependency additions.

Independent NumPy P1 calculations ran twelve allocation cases spanning single
and repeated trajectories, constant/nonlinear materials, matching/nonmatching
contacts, and radiating/nonradiating models, with and without gradient guidance.
A separate Brent root checked every bounded conditional thermal decision. The
largest selected-watt discrepancy was 0.000733 W, below the declared 0.001 W
bracket tolerance, and replay of each selected reference field was bit-identical.
These are independent mathematical references, NOT Rust executions.

For the example, the reference selects standby=3 W, memory=10 W and
chip=14.734143 W. Its all-cycle sampled peak is 302.49999985 K at 20 seconds.
Guidance and bisection each used 19 complete evaluations on this case; the reverse
calculations add cost, so no speedup is claimed. Reversing priorities gives about
16.569090 W to chip and near-zero to the other controls. With only two evaluations,
the reference retains standby=3 W and the other controls at zero, explicitly
unfinished. A common-multiplier search instead keeps standby at zero and preserves
the original chip/memory ratio; it cannot represent this independent allocation.
