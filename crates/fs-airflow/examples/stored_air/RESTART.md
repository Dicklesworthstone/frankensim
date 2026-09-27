# Stop and resume a stored-air FEM run

Build once and retain that executable. A checkpoint binds its exact bytes, the
mesh and vertex ordering, material/source constants, duty history, air capacity,
ventilation, temporal estimator tolerance, and the runtime's complete numerical
policy. A different executable or model refuses rather than resuming a different
problem. The same physical flags and duty file must be supplied on resume.

```bash
cargo build -p fs-airflow --example stored_air

target/debug/examples/stored_air --ventilation-w-k 0 \
  --attempts 4 --checkpoint-dir checkpoints-a > prefix.csv
# A short work budget returns nonzero with an accepted partial trajectory.

target/debug/examples/stored_air --ventilation-w-k 0 \
  --resume checkpoints-a --checkpoint-dir checkpoints-b > suffix.csv
```

`--resume` accepts an individual `.fscp` file or a checkpoint directory. Directory
resume selects the last completed generation, never a `.pending` file. A corrupt
completed generation refuses; it is not silently replaced with an older state.
`--checkpoint-dir` must be a NEW directory, including on resumed invocations.
All previous generations remain untouched. Every completed attempt, including a
rejected attempt's shorter retry duration, is saved before the next attempt.

`--attempts` is an explicitly additional per-invocation allowance. Cumulative
attempt, rejection and producer-evaluation counts survive restart, along with
the entire solid field, air temperature, source energy and net input energy.
A completed checkpoint resumes with no further producer calls. No checkpoint
I/O or executable hashing is performed when neither restart flag is present.

Checkpoint publication writes and syncs a new pending file, then renames it in
the new writer-exclusive directory; Unix also syncs the directory. Durability
still depends on the filesystem. The 1 MiB per-file cap and generation/scan caps
are checked before unbounded reads. A failed write stops the command and leaves
previous completed checkpoints available.

The resumed CSV contains only newly accepted rows and its own header. Keep the
previous CSV for a full trace; checkpoints recover numerical state, not an
atomic external stdout log. For a clean budget stop, concatenate the prefix
with the suffix after removing the suffix header. A hard process interruption
can leave stdout behind the checkpoint and lose the in-flight attempt's work.

The digest checks integrity, not signer authenticity or physical validity.
This remains the spatial FEM / well-mixed-air example, not native `.fsim`/ledger
integration, a CFD model, or a continuous-peak/error certificate.
