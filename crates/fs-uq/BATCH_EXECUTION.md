# Bounded parallel physical uncertainty evaluations

`UqBatchExecution` dispatches at most 256 missing Monte Carlo inputs to a
caller-owned batch driver. The driver may use a TilePool, scoped threads or a
remote evaluation service, but must drain all jobs before returning. It receives
immutable `UqBatchInput` values and one initially empty outcome slot per input.
Fill the matching slot with `Some(Ok(response))` or `Some(Err(diagnostic))`;
leave it `None` only for unlaunched or genuinely unfinished work. Do not suppress
an inconvenient result or error. Catch worker panics as failures, not missing data.

Input-law admission, Philox addressing, singular Gaussian factors and the final
statistics belong to the existing `UqExecution`. The worker count and completion
order never change sample coordinates or reduction order. No additional random
generator, physical model or statistical reduction is introduced.

Each advance is one bounded wave. An unfinished member prevents opening a later
wave, but does not erase successful siblings. The next advance regenerates only
missing inputs at their original ordinals; it can use a different width. Cloning
and the versioned `checkpoint`/`restore` API retain these holes and successful
siblings. Neither restoration nor input regeneration invokes physical callbacks.

`UqBatchReport` separates attempted, accepted, contiguous and out-of-order work.
Only a complete fixed allocation exposes `result`, the original scalar report.
A value-dependent early finish or timeout is not a valid selected-worker sample.
The earliest returned error/nonfinite result permanently refuses the experiment,
even when another sibling succeeded or was interrupted. Refused studies cannot
be checkpointed as resumable prefixes. Callback panics propagate unless the
batch driver contains them; do not resume an execution after an unhandled panic.

Storage is capped at the admitted million observations, 256 pending values and
256 input vectors, each with at most 256 parameters. Pre-dispatch reservation
covers result publication, including after cloning. Physical solver storage and
runtime require separate admission by the driver. Cancellation is checked during
input preparation and after the drained wave; expensive callbacks must implement
their own cancellation. Interrupted costs are driver-owned, as in scalar UQ.

`FSUQBT01` checkpoints bind the complete original plan, sampler semantics and
caller-supplied physical-model identity. Exact response bits, including signed
zero, are retained. The canonical pending window starts at its first missing
ordinal. Bounded framing, statuses, tags, checksum and identity are checked on
restore. Checksums are not authentication: use trusted checkpoint producers and
the same deterministic runtime profile. Scalar checkpoint bytes are unchanged.

`tests/product_batch.rs` exercises real scoped workers, scalar bit parity,
singular/mixed input laws, saved holes, cancellation, failure and malformed files.
Native compilation/execution is required before claiming these tests pass.
No speedup, optional-stopping inference or physical-model validation is implied.
