# Pause and resume a projected multi-load study

The existing `fs-marquee-elasticity-robust --projected` command can now retain
restart checkpoints. This continues the SAME study, unlike the initial-field
CSV option described in `PROJECTED_MULTILOAD.md`, which starts a NEW study.

For a two-case `loads.csv`, start a study and pause after one accepted update:

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-robust -- \
  --projected /tmp/study-part1 loads.csv 4 30 0.45 6 worst 362 \
  --stress-limit 100 --pause-after 1
```

A deliberate pause returns exit **14**, with `status: "paused"`. It is not a
numerical failure, convergence, or completion. The search can instead stall or
exhaust its original budget before the requested pause point. `--pause-after 0`
retains the feasible baseline without attempting an update. `--checkpoint`
alone saves checkpoints without requesting a pause. Neither option changes the
objective, material area, stress policy, or acceptance test.

Continue with a new output directory and an explicit recovery-work allowance:

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-robust -- \
  --projected --resume /tmp/study-part1/checkpoint.fscp /tmp/study-part2 \
  --recovery-solves 4
```

Resume needs **two solves per declared load case**: the original feasible
baseline and the last accepted geometry are independently re-solved. The same
stress sampler rechecks both complete families when a limit is installed. A
bitwise discrepancy with stored numerical evidence refuses recovery. Recovery
solves are reported separately as `checkpoint.recovery_solves_started`; the
original `solves_started` and `max_solves` remain unchanged. Failed recovery
reports its actually started solves on stderr. An exhausted study stays
exhausted; recovery does not replenish its optimization budget.

A resumed invocation may also use `--pause-after N`, counting NEW accepted
updates in that invocation. The original total update goal, global nucleation
schedule, search multiplier, loads and ordering, material parameters, fixed
boundary nodes, area policy, stress limit, candidate limits, original baseline,
and recorded study work all survive. Resume does not accept replacements for
those declarations. The executable retains its grid/update/candidate caps.

## Files and durability

Checkpoints use a bounded versioned binary layout and the existing ledger hash
for byte integrity. `checkpoint-NNNNNN.fscp` is saved at the invocation start and
after each accepted update; `checkpoint.fscp` is saved at a normal numerical,
budget, search, pause or iteration-limit stop. Each file is exclusively created,
flushed and synced. No checkpoint or existing output directory is overwritten.
The summary is written only after final geometry and checkpoint writes succeed.
A write failure does not print success. Earlier completed checkpoint files are
left untouched by later failures; a partially written file fails its hash/length
checks. This is not a crash-atomic transaction over the entire output directory.

Each new output directory contains only its segment's trajectory/attempt rows.
Keep earlier directories for the full trajectory. `baseline-level-set.csv`
remains the ORIGINAL feasible baseline; a resumed `input-level-set.csv` is the
accepted field at this invocation's start. Final metrics and stress remain bound
to the final accepted field, not to the last rejected candidate.

## Limits and verification

Exact continuation requires the same numerical implementation and platform.
A changed implementation must explicitly version or migrate the checkpoint;
there is no cross-version or cross-ISA continuation claim. Checksums do not
prove issuer identity, historical execution, or an anti-replay resource quota.
Abrupt process death can lose work accounting since the last completed snapshot.
Recovery does not serialize tentative candidates; an interrupted candidate is
retried from the accepted boundary with its recorded attempted work preserved.
Recovery solves and existing individual assembly/solver calls are synchronous,
not preemptible or subject to a hard wall-time guarantee.

The new library tests exercise both aggregates, stress/no-stress recovery,
accepted-update split replay, corruption, malformed data and work accounting.
The actual-binary tests cover a stress-constrained split run, exhausted budgets,
invalid checkpoints, immutable policy, and non-overwrite. Focused commands:

```sh
cargo test -p fs-topols --lib robust_descent::engine::projected::checkpoint
cargo test -p fs-marquee --bin fs-marquee-elasticity-robust --test projected_resume --test projected_multiload
```

These Rust tests were authored but not executed in the implementation session:
Cargo, rustc, rustfmt, DSR and RCH were unavailable. Byte-layout reference checks
and source checks are not compilation or a successful elasticity study. This
feature makes no new stress/volume certificate, optimality, physical validation,
or full Journey B completion claim.
