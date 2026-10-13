# Recover a measurement study after process restart

`MeasurementStudy::to_checkpoint(cx)` encodes the original typed alternative
menu, decision threshold, reading cap, ordered actual readings and retained
input/output history. Store the complete returned bytes in durable storage.
`MeasurementStudy::from_checkpoint(&bytes, cx)` reconstructs a fresh study,
recommends each historical request again, applies its actual reading through the
existing Kalman/EVPI path, and checks the retained request and final history.
It never calls a measurement provider and never installs a serialized posterior.

```rust,ignore
let bytes = study.to_checkpoint(cx)?;
// Durably store bytes in the application's file/blob/ledger transaction.
let restored = MeasurementStudy::from_checkpoint(&bytes, cx)?;
let next = restored.recommend(cx)?;
```

All frames have an explicit type/version, canonical little-endian fields, exact
fs-qty quantity-schema tokens and a domain-separated checksum. The decoder
bounds the complete input to 2 MiB before hashing or allocation. It admits at
most 256 alternatives and 4,096 readings, requires sorted unique alternatives,
checks string lengths and finite canonical numbers, rejects trailing/truncated
bytes, and revalidates semantic quantities and relative-versus-physical costs.
The checksum detects accidental corruption; it is NOT authentication. A party
able to replace a complete coherent history can construct another valid study.

Transport, initialization, every replayed recommendation and every assimilation
share one OED work/byte/poll accountant. The existing quadratic action-evaluation
cap applies to all replay rounds together. A valid long history with a large
menu can therefore require more replay work than the synchronous API admits.
Recovery refuses without publishing a partial study, leaving the borrowed
checkpoint available for retry or inspection. There is no unmetered fast path
that trusts stored means or variances.

This is same-build/compatible-numerical-profile recovery. Recomputed request
identities and final history must match exactly; changes in planning, quantity
semantics, Kalman arithmetic or numerical profile may refuse old data. No
migration, cross-ISA guarantee or model/calibration certificate is issued.
Assimilation Colors are recomputed for the recovery execution and may differ
with its budget or context. They are not presented as the original run's
execution evidence; the original readings, numerical state and history remain
bound to the checkpoint. Saving again produces the same scientific transport.

## Keep a completed acquisition across a restart

Save the source study before acquisition. Then durably save
`acquired.to_checkpoint(cx)` BEFORE applying that reading. After restart:

```rust,ignore
let mut restored = MeasurementStudy::from_checkpoint(&study_bytes, cx)?;
let acquired = AcquiredMeasurement::from_checkpoint(&restored, &reading_bytes, cx)?;
acquired.apply(&mut restored, cx)?;
let next_study_bytes = restored.to_checkpoint(cx)?;
```

Recovery checks the source history, current recommended request, exact objective
schema and instrument identity. Neither recovery nor `apply` calls a provider.
Wrong-source and already-applied readings refuse, so an old pending file cannot
silently update the latest generation a second time. The two transports have
different type tags/checksum domains; a study file is not a reading file.
Only valid request-bound acquisitions can be encoded; provider-protocol errors
remain available in their original owned error values for inspection.

An `AcquisitionError::After { acquired, .. }` can be saved using a fresh context
after post-I/O cancellation. Using the cancelled context to serialize still
refuses; no cancellation bypass is added. Keep the original acquired value on
storage failure. Once a complete pending frame has been saved, numerical retry
does not need another sample. A crash BEFORE that durable save is not solved by
this API; device-side idempotency/journaling is needed for that gap.

## Separate-process file example

The example uses a labeled dimensionless prior and caller-supplied outcomes. It
never fabricates a measurement. Each write uses `create_new`, `write_all`, and
`sync_all`, plus parent-directory synchronization on Unix. It never overwrites
or deletes the source study, a reading file, or an existing output generation.
Partial files fail checksum/length admission; retain the last verified source
and retry to a new path. There is no automatic "latest" pointer or multi-file
transaction. Filesystem durability guarantees remain the storage provider's.

```sh
cargo run -p fs-oed-e2e --example recover_sensorforge -- new state-0.fs 8
cargo run -p fs-oed-e2e --example recover_sensorforge -- next state-0.fs
# Supply the alternative printed by `next`, its actual value and instrument:
cargo run -p fs-oed-e2e --example recover_sensorforge -- \
  record state-0.fs reading-0.fs ALTERNATIVE 0.2 actual-meter
# This may run in a completely new process and never invokes the instrument:
cargo run -p fs-oed-e2e --example recover_sensorforge -- \
  apply state-0.fs reading-0.fs state-1.fs
```

The core performs no filesystem I/O and cannot make a caller's write durable.
It can be used with a ledger/blob transaction without an upward dependency from
this numerical crate. These transports do not persist device retry state or
claim exactly-once hardware acquisition, authenticity, or recovery of bytes
that were never successfully stored.

```sh
cargo test -p fs-oed-e2e --test measurement_checkpoint
cargo test -p fs-oed-e2e --test measurement
cargo check -p fs-oed-e2e --example recover_sensorforge
```
