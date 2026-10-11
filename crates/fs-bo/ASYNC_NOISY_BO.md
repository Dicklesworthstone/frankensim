# Keeping external simulation workers busy

`fs_bo::async_study::AsyncNoisyStudy` issues replacement designs without waiting
for a batch's slowest simulation. Construct it with a study ID, dimension,
initial Sobol count, **total job budget**, maximum unanswered jobs, and the
existing `NoisyBoConfig`. The fixed kernel/prior and individual observation
variances retain their original physical units. No new dependencies are used.

```bash
cargo run -p fs-bo --example asynchronous_noisy_design
cargo test -p fs-bo --lib pending
cargo test -p fs-bo --lib async_study
```

`ask()` reads unanswered requests; it does not dispatch or lease them. Keep a
caller-owned map of in-flight job IDs to prevent duplicate external execution.
`tell(request, observation)` records one result. It validates exact request
coordinates, study ID, finite values and nonnegative variance before mutation;
identical retries are idempotent and conflicting responses refuse. Then call
`refill()` to model available results and issue replacements. Initial designs
larger than the concurrency cap are released in order as slots free up. Once
all initial designs have been **issued**, adaptation needs only one received
result; unfinished initial jobs remain pending, not a batch barrier.

Every reservation consumes the total job budget immediately, including jobs
awaiting an external retry. Unanswered reservations count against the in-flight
cap whether dispatched or not. A refill issues at most `config.q` adaptive jobs
and never exceeds either cap. Empty `ask()` is not finalization: after the last
result, `refill()` must still update the final model before `is_complete()` is
true. A repeated refill with no changed data and no capacity does no model work.

## What is optimized

The GP is fitted only to received observations, ordered by logical job ID.
For each joint latent draw, the acquisition computes

```text
max(0, min(latent received designs, latent pending designs)
       - min(latent proposed designs))
```

Thus a pending design's uncertain result is integrated rather than replaced
with its posterior mean or a fabricated observation. This is incremental noisy
expected improvement, following the pending-as-baseline construction described
in BoTorch's incremental qLogNEI documentation. This implementation uses the
existing ordinary q-NEI, not log smoothing, and inherits its numerical jitter,
degenerate-covariance fallback and finite-bank approximation. Exact duplicate
proposals are refused rather than interpreted as valuable replicates. A bounded
search can refuse even when unsearched novel designs exist; it never silently
spends an extra acquisition budget to recover.

The reusable surfaces are `pending::q_pending_noisy_expected_improvement` and
`NoisyBoConfig::propose_pending`. They also support custom external dispatchers.
`AsyncModelReport` names the exact received job IDs used by the last successful
fit and maps its recommendation back to the logical job index. All original
results, including results received since that model update, remain in `jobs()`.
A reported posterior is a model estimate, not a statistical certificate.

## Cancellation and replay

`refill_controlled` polls before and after fitting and at existing greedy-slot
boundaries. A cancellation or model/acquisition refusal does not publish a new
recommendation or partial proposal; issued requests and received results remain
intact. Attempt counters include discarded numerical work. A dense fit or one
slot's CMA-ES search is not preemptible. External solver cancellation, persistent
job leases and failure retries belong to the application; no fake failure value
is inserted into a GP and no failed physical run is counted as a valid result.

Clone provides an in-memory checkpoint. Replay requires identical job outcomes
**and identical availability sets at decision time**. Reordering `tell` calls
without a refill between them does not change a decision, but running a refill
after different subsets arrive can intentionally change the trajectory. The
example illustrates unequal durations with a deterministic logical-time event
queue, not sleeps or real workers. Its first job is a straggler while other
slots continue through an 18-job budget. This is not a wall-clock benchmark.

This asynchronous API is fixed-kernel and unconstrained. Learned/constrained
asynchronous models, durable asynchronous checkpoints, authenticated dispatch,
delay/censoring models and value-of-information replicate selection are not
added. Existing fixed/learned batch APIs and their checkpoint formats are
unchanged. Model uncertainty is not a claim of optional-stopping validity when
runtime or missingness depends on objective values.

Four acquisition regressions cover pending-free parity, joint-draw marginal
identities, near-pending covariance, duplicates, bounded selection and replay.
Nine driver regressions cover stragglers, logical-ID posterior mapping,
initial/final caps, availability-order replay, cancellation, response retries,
model failures, one-worker synchronous parity and pre-work admission. Native
compilation, these tests and the example remain unverified in the authoring
runtime, where DSR/RCH/Cargo are absent. Independent Gaussian equation and
logical-scheduler checks are not execution of the Rust implementation.
