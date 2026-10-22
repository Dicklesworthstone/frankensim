# Concurrent physical cooling uncertainty

`cooling-network-uq BASE UQ --workers N` runs up to N simultaneous original
cooling child processes, where N is explicit and in 1..64. No alternate thermal
solver, approximate response or random input generator is introduced.

```sh
frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json --workers 4
```

The sample budget counts physical evaluations, not waves. Every input has the
original scalar Philox seed/ordinal; successful responses enter the original
statistics in ordinal order. For a deterministic physical callback on the same
runtime profile, the complete result has the same schema and arithmetic as the
serial command. Worker count and completion order do not enter the random law.
Steady and complete-transient sampled-peak observations are both supported.
One transient sample remains one full trajectory, not a timestep or partial peak.

## Keep completed work across interruption

```sh
frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json --workers 3 \
  --checkpoint first.batch --max-new-samples 3

frankensim --json cooling-network-uq \
  examples/cooling-network/fan-correlated-hotspot.json \
  examples/cooling-network/uq-fan-hotspot.json --workers 1 \
  --resume first.batch --checkpoint complete.batch
```

The first invocation exits BUDGET with progress only. The second uses a different
worker count and evaluates only missing samples under the original lifetime cap.
`--workers 1` is the serial execution mode for this checkpoint format. Ordinary
scalar checkpoints and `FSUQBT01` batch checkpoints are deliberately distinct;
repeat --workers when restoring a batch study. Model, input laws, seeds, limits,
executable and total budget stay bound. Worker count, invocation wall allowance
and chunk size may change. Checkpoints must be trusted, not merely checksummed.

A completed sibling behind an interrupted earlier sample is retained in its
original slot, including on disk. No later wave opens until those gaps finish.
Progress reports distinguish `contiguous_observations` and
`retained_out_of_order`; neither subset yields statistics. No failed response
is skipped. Genuine physical errors, nonfinite outputs, failed worker creation
and unwinding worker panics refuse the experiment after draining launched jobs.
The current checkpoint output becomes a terminal failure diagnosis. Inputs are
never overwritten and each output path must be new.

Each invocation uses one wall deadline shared by its children. The existing
watchdog kills and reaps interrupted children; scoped producer threads are joined
before a wave is recorded. All returned complete values survive a later timeout.
A completed checkpoint returns its result without new physics. Explicit zero
chunks launch nothing. Checkpoint publication happens after each drained wave,
not individually mid-wave: a process crash can lose the current unsaved wave,
but cannot destroy the previous atomically published checkpoint.

## Resource and inference scope

This mode currently admits fixed-count Monte Carlo only. QMC, importance, PCE,
sensitivity, mean controls, sequential decisions and candidate selection conflict
with --workers and refuse before input reads. Existing modes without --workers
are unchanged. The CLI still admits at most 10,000 physical evaluations.

The process cap is not a global CPU/RAM limit: each child can itself use solver
threads and memory, and the existing bounded stdout/stderr buffers are per child.
Choose N to fit aggregate resources. No measured speedup or memory bound on the
physical solver is claimed. No partial distribution, confidence sequence,
continuous-time peak certificate or physical qualification is added.

Focused native regressions:

```sh
cargo test -p fs-uq --test product_batch
cargo test -p fs-cli --bin frankensim uq_command::execute::parallel
cargo test -p fs-cli --test cooling_parallel
```

These tests must run natively before claiming Rust execution or bit-parity passes.
