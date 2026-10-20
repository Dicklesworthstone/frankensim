# Long scan objectives with short adjoint tapes

The full-storage scan objective keeps one sparse thermal linearization per
endpoint. Checkpointed evaluation instead retains caloric states at segment
boundaries, replays one segment at a time through the same physical solvers,
and drops that segment's tangents after its reverse sweeps finish.

This is uniform segmentation, not an optimal Revolve scheduler. It trades an
additional forward pass for a bounded number of simultaneously live sparse
tangents. It does not replace the heat equation, radiation law, source
projection, phase chart, nonlinear solve or physical acceptance tests.

## Call the physical evaluator

The module is
`fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::checkpoint`.

- `evaluate_scan_objective_checkpointed` returns the original
  `ScanObjectiveEvaluation` plus `ScanCheckpointWork`.
- `evaluate_scan_objective_checkpointed_with_peaks` returns that same objective,
  independent point responses, sampled hot-spot envelopes and replay work.
  Each segment is replayed once for all outputs, not once per constraint.

Both functions consume the trial `ScanningEnthalpyRun`, just as the existing
full-storage objective does. Every checkpoint retains specific enthalpy,
temperature, liquid mass fraction, time and logical step count together. The
mesh, reference masses, phase charts, fixed contacts, conductivity assignments,
source schedules, radiation and motion policy remain in the same run. There is
no public rewind operation on application-owned state.

For N endpoints, V vertices, S sources and effective segment size
B = min(segment_steps, N), admission requires:

| Allowance | Meaning in checkpoint mode |
| --- | --- |
| `ScanObjectiveBudget::max_steps` | At least N, still the whole trajectory |
| `max_retained_vertices` | At least B * V |
| `max_retained_source_vertices` | At least B * S * V |
| `max_retained_radiation_entries` | Combined feedback allowance for one segment |
| `ScanCheckpointBudget::max_checkpoint_vertices` | At least (ceil(N / B) + 1) * V |
| `ScanCheckpointBudget::max_forward_solves` | At least 2 * N |

Each checkpoint vertex represents three scalar arrays, not one scalar or one
byte. Existing reports, sparse matrices, solver workspaces and allocation
overhead remain additional memory. Peak responses retain their selected scalar
samples under `max_peak_samples`; their running gradients use the existing
response-gradient allowance. These are shape/work caps, not peak-byte claims.

For example, N=4096 and B=64 means at most 64 endpoint linearizations coexist,
with 65 caloric checkpoint states, rather than retaining 4096 linearizations.
This is a storage-shape illustration, not a measured process-memory reduction.
B=1 minimizes tangent storage but keeps N+1 caloric snapshots; it is not
necessarily the smallest overall memory choice.

Successful evaluation reports N original and N recomputed primal endpoints.
Adding point or peak constraints adds their admitted transpose work, but no
additional forward replay. A failed evaluation returns the original typed
error, not a partial objective or gradient. Newton/Krylov work remains bounded
by the original per-endpoint policies; endpoint counts are not iteration counts.

## Use it with optimization

Call `ScanPowerStudy::set_checkpointing(cx, policy)` before any physical probe,
including optimizer initialization or construction of `ConstrainedScanStudy`:

```rust
use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::checkpoint::ScanCheckpointBudget;

study.set_checkpointing(cx, ScanCheckpointBudget {
    segment_steps: 64,
    max_checkpoint_vertices: 65 * vertices,
    max_forward_solves: 8192,
})?;
```

The example numbers above describe 4096 endpoints. The original objective's
retained-tangent allowances must separately cover the selected segment. For
runtime-computed sizes, use checked integer arithmetic before constructing the
policy. Mesh-specific admission occurs before the evaluator runs physics.

The selected policy applies to L-BFGS, SQP, objective replay, point-state
constraints and hot-spot constraints. A study cannot change it after a physical
evaluation has been reserved, including a failed probe. Extending an optimizer's
evaluation budget neither changes the memory policy nor refunds earlier work.
`ScanPowerEvaluation::checkpoint_work` reports the last evaluation's actual
successful original/recomputed endpoint counts and retained shapes. `None`
means the original full-storage path was used. Source export still spends no
physical evaluation.

The standalone synthetic fitting example accepts an even endpoint count and a
segment length:

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --checkpointed
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --checkpointed 64 8
cargo test -p fs-conduction --test scan_checkpoint
cargo test --manifest-path examples/scan_power_fit/Cargo.toml
```

The default checkpointed example uses 64 endpoints and segments of 8, fitting
the same two held heating levels as the original example. Its final replay
prints 64 original endpoints, 64 recomputed endpoints, 9 checkpoint states and
8 peak live tangents. Those diagnostics refer to the final replay, not the
cumulative optimization workload. Unconverged fits still exit unsuccessfully.

## Replay and numerical scope

Every replayed segment endpoint must match all retained caloric arrays, time
and logical count bit-for-bit. Replayed fitting observations, point responses
and selected peak samples are also compared bitwise before publication. Full
storage and replay use the same reverse accumulator and peak reduction, in the
same chronological/summation order. A mismatch is a refusal, not a reason to
relax solver tolerances. This checks the reached replay states and selected
samples; it is not exhaustive proof of every execution interleaving.

Time grids and event sides remain fixed. Event-left grids must include all
source, boundary, reservoir and selected motion stops; omitted stops refuse.
Replay never silently changes the grid. Material/phase corners still refuse
classical derivatives. Peak envelopes retain their explicit smoothing margin
and concern selected nodes/times only, not unsampled physical maxima.

This is in-process checkpointed differentiation, not durable serialization,
remeshing, adaptive-grid differentiation, a global error bound or a safety
certificate. The authored Rust regressions exercise full-storage equivalence,
short-tape admission, held events, existing history, cancellation, combined
point/peak responses and both optimizers. Their native compilation and execution
were not established in the authoring environment because Cargo/rustc and the
configured DSR/RCH runners were absent. Independent NumPy/SciPy equation checks
are separate evidence, not execution of the Rust implementation.
