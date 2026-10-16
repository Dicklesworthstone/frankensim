# Pausing and resuming a stress-constrained study

The single-load normalized 2-D cantilever supports durable continuation through
`fs-marquee-elasticity-stress --projected`. The original whole-trajectory mode
and ordinary projected invocation remain available.

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-stress -- \
  --projected ./stress-part-1 1e12 3 30 0.6 8 0 300 --pause-after 1
# Exit 6 means the requested pause, not a completed/converged study.

cargo run --release -p fs-marquee --bin fs-marquee-elasticity-stress -- \
  --projected --resume ./stress-part-1/checkpoint-0001.fscp \
  ./stress-part-2 --wall-seconds 300
```

The deliberately loose `1e12` sampled-stress limit above is a software example,
not an engineering design limit. `--checkpoint` enables checkpoint retention
without requesting a pause. `--pause-after N` implies checkpointing and counts
accepted updates in the **current invocation**; `0` retains the feasible baseline.
Resume also accepts `--pause-after N`, but cannot override any physical, stress,
volume, candidate-search, or total-iteration setting.

Each checkpoint retains exact floating-point node bits, the global update
ordinal, augmented-Lagrange multiplier, complete optimizer controls, fixed
support/load traces, area projection policy and sampled-stress limit. A complete
checkpoint is flushed before more physics; incomplete writes fail digest or
length admission and do not replace earlier checkpoints. Resume writes a NEW
output directory and never changes the source study. Retain the original binary:
a different executable is refused, even if its package version is unchanged.
The bounded reader admits at most 1 MiB, levels 2..=7, and 200 total updates.

Recovery independently solves compliance/area and sampled stress before it
publishes a new study. Saved metrics must agree bit-for-bit. The digest is a
corruption check, not authentication or proof of historical execution. Replay
assumes the same admitted deterministic execution profile; no cross-ISA promise
is added.

The wall budget covers baseline construction, checkpoint recovery and subsequent
updates. It is polled inside the existing true-residual CG solves, between area
projection evaluations, and before each cell's stress probes. The retained
`poll_iters` controls CG batching in both proposal and final solves. Assembly,
individual area-quadrature evaluations, one cell's probes, bounded checkpoint
I/O/decoding and executable fingerprinting remain indivisible. This is a
cooperative work boundary, not a hard wall-time deadline.

A timeout before full baseline/recovery admission, or before the first study
publication, exits 6, writes a diagnostic to stderr and creates no output study
or success summary. The source checkpoint is untouched. A timeout after study
publication retains the last fully accepted field/checkpoint and reports
`wall_budget`; it never publishes the interrupted candidate or a partial stress
maximum. Retrying recovery re-solves the saved state, not its previous updates.

Summary schema `fs-marquee-projected-stress-v2` names the last checkpoint,
global `accepted_updates`, `start_iteration`, and `segment_accepted_updates`.
The baseline and reduction are explicitly scoped to the current segment;
previous history remains in the source study. `iteration_limit` exits 0,
`paused`/`wall_budget` exit 6, and `no_feasible_descent` exits 11. None implies
convergence, KKT optimality, a continuous stress bound, or physical validation.
