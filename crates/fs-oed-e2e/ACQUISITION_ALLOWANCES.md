# Bound acquisition spending before contacting an instrument

`plan::acquisition::budget::AcquisitionAccount` holds a finite cumulative
allowance and named-attempt cap. It is separate from the `Cx` numerical-work
budget and the study's accepted-reading cap. A failed physical attempt may
consume acquisition resources without contributing an accepted reading.

```rust,ignore
use fs_oed_e2e::plan::acquisition::budget::AcquisitionAccount;

let mut account = AcquisitionAccount::new(total_allowance, max_attempts)?;
let permit = account.reserve(&study, &request, "experiment-42-attempt-0", cx)?;
// Durably save the source study and this debit BEFORE contacting the instrument.
let account_bytes = account.to_checkpoint(cx)?;
// Save account_bytes in the application's durable store here.
let acquired = permit.acquire(&study, &mut instrument, cx)?;
// Save acquired.to_checkpoint(cx) BEFORE numerical application.
acquired.apply(&mut study, cx)?;
```

The allowance and sensor costs must have the exact same `AcquisitionCostSpec`,
including the relative-weight versus physical-quantity role. Choose an additive
resource such as seconds or joules; relative weights remain abstract points,
not measured money or energy. The amount is finite and positive; setting the
attempt cap to zero disables acquisition. At most 4,096 attempts are retained.

Reservation validates the current study request, attempt identity, cost schema,
count and cumulative amount before changing the account. Failed admission does
not call a provider, change the study, or consume allowance. Each successful
reservation returns a non-Clone permit consumed by at most one provider call.
The ordinary acquisition path rechecks staleness and cancellation before I/O.

A commitment is NEVER automatically refunded: dropping a permit, provider
failure, stale-request refusal, post-I/O cancellation and numerical failure do
not restore the allowance. `charges()` records commitments, NOT successful
measurements or invoices. The requested cost is a caller-declared resource
charge, not a measured bound on arbitrary provider behavior.

An attempt ID cannot be reserved again. A deliberate retry must name a new
attempt and consumes another charge. This avoids treating a repeated application
request as free hardware work. Applying an existing `AcquiredMeasurement`
(including one retained in `AcquisitionError::After`) neither calls a provider
nor reserves another attempt.

## Restart without resetting the allowance

`AcquisitionAccount::from_checkpoint(original_allowance, original_attempt_cap,
&bytes, cx)` verifies the configuration and replays the bounded charge history.
It reconstructs the total rather than trusting a saved balance, rejects repeated
attempt IDs and costs outside the allowance, and performs no I/O. Save the account
after reservation, before using its permit. Permits cannot be restored: when I/O
status is unknown after a crash, recover a saved reading or explicitly pay for a
new attempt. There is no inferred cancellation refund or implicit re-acquisition.

The account frame is independent of the scientific study checkpoint. Both are
needed after restart: accepted readings alone cannot reconstruct failed or
abandoned acquisition spending. The caller must choose the latest durable
account generation and coordinate writers. A checksum detects corruption; it
does not authenticate a ledger, prevent deliberate rollback/copying, or provide
a distributed spending lock. The library cannot recover a debit never saved.

## Numerical and planning scope

Cost accumulation rounds conservatively upward only when ordinary binary64
addition rounds below the exact sum. A positive tiny cost cannot disappear in a
large accumulated balance. `committed_cost_upper()` is an upper bound on declared
commitments and can conservatively refuse near the cap; it is not an exact bill.
This enforces acquisition limits without changing the scientific recommendation.
An unaffordable winning request is refused, not silently replaced by another
sensor, and exhausting the allowance does not establish decision robustness.

Each reservation or checkpoint operation has bounded allocation and cancellable
work. No runtime dependency or I/O backend is added.

```sh
cargo test -p fs-oed-e2e --test acquisition_budget --test measurement_checkpoint
cargo test -p fs-oed-e2e --test measurement
```
