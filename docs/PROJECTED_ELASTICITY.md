# Hard-area-constrained elasticity study

The `--projected` mode of `fs-marquee-elasticity` optimizes the existing normalized
unit-square, plane-strain cantilever using a free-boundary bilinear level set.
It reuses `fs-topols` evolution and the canonical CutFEM elasticity solver. The
existing unflagged whole-trajectory guarded mode is unchanged.

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity -- \
  --projected /tmp/projected-study 4 30 0.45 6
```

The positional arguments are output directory, grid level, requested accepted
updates, material-area target, and candidate budget per update. Defaults are
level 4, 30 updates, area 0.45 and six candidates. The output directory must not
already exist. The driver bounds levels to 2–7, updates to 1–200, and candidates
to 1–16. The area tolerance is `1e-4` of the normalized design-box area. Projection
uses at most 64 area evaluations and an absolute free-node offset of at most 2.

## What changes

The supplied strip first undergoes **feasibility restoration**. Its projected,
independently solved geometry becomes the objective baseline. Material removal
from an overfilled strip is not presented as a compliance improvement.

Every update starts from the current feasible geometry. The existing engine
proposes an evolved field; the new projection imposes the area equality; the
canonical elasticity operator then solves the **projected** field. Acceptance
requires both the declared area tolerance and a strict relative compliance
reduction (default `1e-8`). Rejected candidates contract the proposal travel by
one half, always starting again from the last accepted geometry. No accepted
update is reported from an unprojected proposal's objective.

Both left and right boundary nodal traces remain fixed at their declared initial
values. This preserves the clamp and loaded boundary geometry instead of making
the objective smaller by removing its load support. Library consumers may
prescribe any strictly ordered set of fixed nodal values through `volume` and
`projected`. A fixed node does not silently freeze the whole adjacent cell.

## Outputs and stop meanings

`baseline-level-set.csv` and `level-set.csv` contain the actual feasible baseline
and final accepted nodal fields. `trajectory.jsonl` contains independently solved
accepted-state objectives, areas and snapshots. `attempts.jsonl` retains failed
area/solve/decrease gates. `unprojected-proposals.jsonl` is separate evolution
information; its audits and objectives do not certify the projected geometry.

`summary.json` is written last, after the fields and traces have been flushed.
A file-write error returns failure and does not print a success result. Accepted
trace rows are flushed as updates complete, not accumulated until the end.

`iteration_limit` (exit 0) means only that the requested accepted-update count
was reached. `no_descent` (exit 11) means the bounded candidate family found no
acceptable next state, including numerical refusals. Neither means convergence
or global optimality. A failed or stalled search leaves the last accepted field
intact. The report retains the measured reduction from the feasible baseline;
there is no asserted 30-percent improvement target or tolerance tuned to it.

## Library continuation and cancellation

`ProjectedOptimizer::checkpoint()` exposes the accepted geometry, global ordinal
and proposal multiplier through the existing `OptimizeCheckpoint` type.
`from_checkpoint` re-solves and admits an already feasible checkpoint without
quietly changing its field; its baseline is explicitly the start of that resumed
segment. The caller retains previous segment evidence. This is a library API,
not a new CLI checkpoint-file format.

`advance_one_controlled` polls the existing proposal CG solves, each complete
area evaluation, and final publication. Interruption never replaces the accepted
state. Final independent elasticity evaluation and individual area evaluations
are checked at boundaries, not internally preempted; there is no hard wall-time
latency guarantee. The proposal's AL multiplier remains heuristic search state,
not a KKT multiplier for the projected constrained problem.

## Verification and scope

The new Rust tests cover analytic area projection, fixed boundaries, refusal and
cancellation, feasible baselines, independently solved accepted endpoints,
continuation, and the real binary's geometry/trace outputs. Suggested focused
commands from a configured constellation checkout:

```sh
cargo test -p fs-topols --lib volume
cargo test -p fs-topols --lib projected
cargo test -p fs-marquee --bin fs-marquee-elasticity --test projected
```

Native execution of these tests is pending: the implementation session had no
Cargo, rustc, DSR or RCH. Independent Python exact-cell projection controls and
lexical/delimiter checks ran, but do not constitute a Rust build or PDE result.
The material area is numerical cut quadrature, not a guaranteed continuum bound.
This is not a 3-D solver, an experimental validation, a certified optimum, or a
completed full-scope Journey B milestone.
