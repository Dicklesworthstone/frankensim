# Discrete transient cooling adjoints

The `cooling-network` command differentiates the selected final temperature
or sampled trajectory peak through its actual backward-Euler solid and mixed-air
endpoint equations. This is separate from the steady gradient fields.

```bash
frankensim --json cooling-network examples/cooling-network/adjoint-contact-pulse.json
frankensim --json cooling-network examples/cooling-network/adjoint-repeated-contact-pulse.json
```

Keep `objective.gradient` false: that field requests a steady adjoint. Add an
explicit policy inside the existing `transient` object instead:

```json
"adjoint": {
  "qoi": "sampled-peak",
  "max_checkpoint_bytes": 1048576
}
```

`qoi` is either `final` or `sampled-peak`. The former uses the selected spatial
objective at the final endpoint. The latter uses its largest value over the
initial state and every accepted endpoint, selecting the earliest exact time
tie. A spatial maximum keeps the existing lowest-ID active-vertex rule. The
result is the derivative of that selected branch, not a unique smooth maximum
derivative at a tie or a continuous-time peak bound. Area-mean objectives also
work, with named regions bound in network order rather than declaration order.

## Repeated duty cycles

A fixed repeated schedule can use the same adjoint policy:

```json
"repeat": {"cycles": 3, "max_total_steps": 21},
"adjoint": {"qoi": "sampled-peak", "max_checkpoint_bytes": 1048576}
```

There are seven accepted steps per cycle in `adjoint-repeated-contact-pulse.json`.
Its complete trajectory has 21 steps and spans 42 seconds. The forward calculation
still uses local cycle time; global reporting offsets do not alter any timestep.
Every cycle starts from the previous accepted field, not from a cold reset.

The adjoint is returned in `repeated_cycles.adjoint`, with global `time_s` and
`state_index`, rather than in `transient.adjoint`. The latter remains null
because `transient` describes only the final cycle in local time. `final`
selects the end of the last cycle; `sampled-peak` selects the largest value
across the ENTIRE trajectory, including warm-up and cycle-start samples.

One reverse traversal propagates through all earlier cycles. It never groups
all occurrences of an interval before reversing the intervening intervals.
Equal interval labels in different cycles identify a shared control, not an
excuse to change chronological order or reset the history multiplier.

The example's independent dense NumPy FEM reference has a sampled peak near
303.8150504 K at global time 34 seconds. The cooldown fan control has a nonzero
derivative at that peak: although the final cooldown has not happened yet,
earlier occurrences of the same cooldown control changed the warm-up history.
These are illustrative numerical model inputs, not measured or validated data.

## Returned controls

For a single cycle, `transient.adjoint` reports the selected value, time, state
index and active vertex. For repeated runs those fields are under
`repeated_cycles.adjoint`. Both report:

- `dtemperature_dinitial_temperatures`: one partial derivative per ORIGINAL
  initial solid node. `dtemperature_duniform_initial_k` sums them for a uniform
  initial shift. The initial field is not reset or independently varied later.
- `dtemperature_dcapacity_multiplier_k`: derivative at multiplier one when ALL
  declared element heat capacities are scaled together throughout the run.
- `dtemperature_dinlet_temperatures`: partial derivatives per hydraulic node for
  supply-temperature changes maintained throughout the run; nonsupplies stay zero.
- Per-interval `dtemperature_dpower_multiplier_k` and
  `dtemperature_dlog_fan_speed_ratio_k`. The first scales that interval's ENTIRE
  currently declared source at multiplier one, with fixed footprints. It is not
  an absolute-watt or individual-component derivative. The second includes the
  existing single-bank affinity, air-capacity and supported convection response;
  it is null without a fan or when the required convection slope is unavailable.

For repeated runs, one interval control changes EVERY occurrence of that base
interval. Its output is the sum of the causal contributions from those
occurrences. An occurrence strictly after the selected peak contributes zero,
but earlier occurrences can still contribute. A zero-source interval has zero
power-multiplier sensitivity. An initial-state maximum needs no reverse
endpoint solves. None of the steady gradient fields is repurposed.

## Adjoint-guided transient sizing

Combine the explicit `adjoint` policy with either `transient.power_design` or
`transient.fan_speed_design` to guide the existing evaluated-bracket search:

```bash
frankensim --json cooling-network examples/cooling-network/size-transient-adjoint-power.json
frankensim --json cooling-network examples/cooling-network/size-transient-adjoint-fan.json
```

Both examples use three fixed cycles of 120-second heating and 120-second
cooldown, with 20-second steps: 36 accepted endpoints spanning 720 seconds per
candidate. The first scales named chip watts toward a 314 K sampled-peak target;
the second scales the complete fan schedule toward 316.4 K. All material, fan,
resistance and power inputs are illustrative declarations, not measured hardware.

Sizing requires `adjoint.qoi = "sampled-peak"`; `final` is explicitly refused
because a final-temperature derivative cannot represent the peak constraint.
The complete all-cycle peak decides feasibility. Each candidate starts at the
original initial field, and its actual speeds and loads are bound into a fresh
schedule used by BOTH forward and reverse solves. No prior candidate field or
adjoint history is recycled. The adjoint summary stays typed through the search;
no rendered JSON is parsed to recover a derivative.

The search sums the shared interval controls once and divides by the candidate's
outer multiplier. Its `history[].dpeak_dmultiplier_k` therefore differentiates
the outer design coordinate at the evaluated candidate. The per-interval
adjoint fields still describe relative changes to that candidate's actual
loads/speeds. At zero workload the relative-to-outer conversion is unavailable,
not a division by zero. Unavailable, wrong-sign, zero or unrepresentable slopes
fall back to bisection. A requested adjoint's solver/domain failure remains an
error, not an infeasible-temperature observation.

Newton proposals must lie strictly inside the central 80% of the already
evaluated bracket. Every proposal reruns the complete physical trajectory and
requested adjoint; predicted temperatures never establish feasibility. The
existing passing-temperature slack and multiplier-width gates both remain
mandatory. `search_method` identifies the policy and `newton_trials` counts
completed derivative-proposed candidates. Without `adjoint`, the search keeps
bisection and performs no hidden reverse solves.

This is opt-in, not a universal speedup: adjoints cost additional work and a
safeguarded search can take more evaluations than bisection. The original wall
budget covers the entire search, including all reverse solves. Existing
per-candidate timestep limits, per-endpoint derivative/linear limits and the
whole-trajectory checkpoint cap are not enlarged. Reconstruction work is in
each trial's `solid_solves` and the design's aggregate `total_solid_solves`.

## Numerical implementation and resources

The endpoint equation is `C(T_new-T_old)/dt + A(T_new)T_new - b = 0`. Its solid
Jacobian is `C/dt + J(T_new)`, including the existing material K-prime and
matching-contact terms. The shared Robin-response algebra and mixed-network
implicit adjoint close the air feedback before the total nodal-load multiplier
is carried backwards as `(C/dt)^T lambda`. Internal Newton/IQN trials never
advance physical history and are not differentiated.

Only accepted temperature and reference vectors are retained. After the COMPLETE
trajectory and its per-cycle/cumulative energy checks pass, the reverse
traversal reconstructs one endpoint at a time using its retained previous
field, local timestep and references. Exact same-profile temperature replay is
required before using that operator. It does not rerun a complete trajectory
per control or retain every solver matrix.

`max_checkpoint_bytes` bounds retained field/reference and frame storage, up
to 512 MiB. The charge covers ALL requested cycles before the first PDE solve,
not a fresh allowance per cycle. It is not a total solver-workspace or allocator
memory guarantee. These checkpoints are in-memory derivative state, not durable
resume files. The original wall deadline covers forward and reverse work;
derivative/linear iteration budgets apply at each reverse endpoint.

Failures or cancellation publish no partial adjoint. Reconstruction solves are
included in the complete trajectory's `total_solid_solves`, with
`forward_solid_solves` separately reported. The adjoint reports
`reconstructed_solid_endpoints`, `adjoint_sweeps` and maximum interface residual.
For repetitions, each cycle summary and the final-cycle `transient` record still
count their own forward work only; complete reverse work belongs to
`repeated_cycles`. Forward nonlinear statistics do not include reconstruction.

Adjoints require a fixed accepted time grid and fixed cycle count. Adaptive
timesteps, `repeat.until_periodic` and hysteretic `repeat.fan_controller` with
adjoints refuse explicitly, including inside sizing. Ordinary non-adjoint runs
keep those existing capabilities. Geometry, material laws, contact resistance
and fluid properties remain fixed. Material derivative kinks or sampled validity
endpoints refuse rather than invent a unique slope. No infinite-cycle,
phase-change, continuous-peak, physical-validation or uncertainty certificate
is inferred. A passing bracket endpoint is not a global optimum or a hardware
compliance certificate.

## Focused checks

```bash
cargo test -p fs-conduction --test backward_euler_adjoint
cargo test -p fs-cli --test cooling_transient_adjoint
cargo test -p fs-cli --test cooling_transient_mean_adjoint
cargo test -p fs-cli --test cooling_repeated_adjoint
cargo test -p fs-cli --test cooling_transient_adjoint_design
```

The repeated regression tests cover complete nonlinear contact trajectories,
finite differences, explicit unrolling and shared-control sums, exact forward
field/history retention, single-cycle equivalence, initial all-cycle peaks,
whole-trajectory checkpoint admission, and unsupported/interrupted runs.
The design tests add actual-candidate replay, non-unit multiplier chain rules,
fan/workload searches, named watts, zero workload, derivative-free controls and
refusals. The authoring environment has no Rust toolchain: these test sources
were not executed there. Independent dense NumPy FEM checks are numerical
references, not Rust execution or runtime benchmarks.
