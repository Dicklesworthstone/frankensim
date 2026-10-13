# Streaming calibration-aware assimilation and sensor design

Use `groups::stream::SharedNoiseStream` when batches share the SAME fixed
calibration or registration error. Calling `JacobianBatch::assimilate` repeatedly
on only the last physical marginal resets independence between the physical
state and nuisance; it is not the correct update for shared errors across calls.

The stream retains `[physical state, nuisance parameters]` in a UD factor,
including learned cross-covariance. It consumes checked `JacobianBatch` values
but applies `[H, J]` with **independent instrument noise only**. Initial physical
state and zero-mean nuisance are independent. This is a static additive Gaussian
model; there is no time propagation, nuisance drift, nonlinear relinearization,
or certification that the supplied calibration/units are scientifically valid.

## Retained updates

Construct `SharedNoiseStream::new(&prior, &model, row_capacity, cx)`, then call
`stream.assimilate(&batch, cx)` for each arrival. `joint()` exposes the complete
read-only factor; `physical_belief(cx)` optionally reconstructs a checked dense
physical marginal WITHOUT discarding the joint state retained by the stream.

A row ID denotes one acquired record, not a sensor type. A genuine replicate
needs a distinct row ID and independently justified instrument noise. Exact
repeats are ignored; changing the value, variance, instrument, operator or
Jacobian under an accepted ID refuses the entire arrival. Updates report new
and replayed batch indices. Numerical diagnostics describe new augmented rows,
not physical-only marginalized residuals. Failures leave the posterior and row
history unchanged. The input-history root is insensitive to chunk boundaries
and exact repeats, but preserves the order of newly acquired rows.

Limits: at most 256 augmented components, 64 nuisance parameters, 256 readings
per checked batch and 4,096 retained row identities, with a caller-selected row
cap. Replay protection is never silently evicted. It is in-memory only; creating
a new stream does not recover a previous history. The existing batch constructor
still performs bounded dense covariance admission. Factor updates cost
O(m (n+k)^2); optional dense publication costs O((n+k)^3). Preparation and the
factor sequence retain separate budget scopes, not one global stream-lifetime
quota. Neither a history hash nor scalar diagnostics authenticate a calibration.

## Choose the next measurement without reading its outcome

`groups::stream::design::SharedNoiseProbe::new` takes a model, fresh row ID,
physical operator, signed nuisance Jacobian and independent variance. It has no
measurement-value or oracle-truth argument. `stream.rank_probes(&target, &probes,
cx)` ranks expected reduction of the variance of the linear physical quantity
`target * x`, using the current joint posterior and actual declared noise.

Every score is one-step and refers to the SAME retained history. Re-rank after
acquiring data; this is not joint batch selection, EVPI, or a cost-aware optimum.
Ties use row ID. A repeated same-direction reading can carry little information
near a shared-error floor, while an opposite-loading measurement can resolve
that ambiguity. The score exposes predicted measurement mean/variance and
before/after target variance. The latter uses a Joseph sum of squares rather
than cancellation-prone prior-minus-improvement subtraction. Unresolved range
or contraction refuses rather than clamping; these are numerical predictions,
not covariance-model certificates. Planning never changes the stream.

After a real acquisition, `probe.record(&model, measured_value, instrument)`
creates the matching checked record for a `JacobianBatch`. Prospective row IDs
already accepted by the stream are refused. A known-zero target variance earns
zero information, not infinite value. Operators are nonzero physical rows, as
required by the existing observation API; calibration-only rows are not added.

Run the synthetic example and all related regression targets:

```sh
cargo run -p fs-assimilate --example shared_stream
cargo test -p fs-assimilate --test factor --test robust --test factor_batch --test shared_noise
```
