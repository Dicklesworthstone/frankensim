# Full shared calibration and registration covariance

`groups::jacobian` extends unit-loading common-mode groups to a named,
multivariate nuisance model. Given independent instrument noise variances
`r_i`, nuisance covariance `C`, and signed observation Jacobian rows `J_i`, it
constructs `R = diag(r_i) + J C J^T` and assimilates through the existing checked
correlated-factor path. No dense inverse or per-sensor variance approximation
replaces the full covariance.

For example, a scan's shared offset and tilt have Jacobian `[1, position]`.
Opposite signs of position matter. A positive common offset limits the
information gained by repetition; opposite signed sensitivities can instead
cancel a nuisance. Correlation between the offset and tilt parameters remains
in `C` and is propagated into every affected pair of observations.

## API

Construct a `SharedNoiseModel::new(identity, parameter_names, covariance, cx)`.
Then bind each `JacobianRecord::new(observation, row_id, &model, jacobian)` to that
model's parameter ordering. `observation.noise_var()` must contain **independent
instrument noise only**. `JacobianBatch::new(model, rows, cx)` checks row
uniqueness, rejects mixed models, and materializes the full observation
covariance. `batch.assimilate(&prior, cx)` returns `JacobianAssimilation`, whose
`robust()` exposes the posterior, numerical audit, and factor-row diagnostics.

The public model, batch, and result content identities bind the ordered
parameter names, full covariance, original instrument/row identities, actual
readings, observation operators, and signed Jacobians to the numerical audit.
Even two Jacobians producing the same covariance retain different identities.
Duplicate prevention is within one batch, not a persistent replay database.

Projection uses `T = J U sqrt(D)` from the nuisance UD factor, then
`R = diag(r_i) + T T^T`. Its cost is O(m k² + m² k), for at most 256 observations
and 64 nuisance parameters. Allocation, projection, identity hashing, and row
checks are bounded/cancellable. This projection has one cumulative accountant;
model validation and robust assimilation retain their checked stage scopes.

Singular **nuisance** covariance is supported when the existing PSD gate and UD
factorization admit it. Singular **observation-noise** covariance is not: the
strict Cholesky gate is unchanged. Numeric overflow or unresolved admission
refuses rather than adding jitter, dropping cross-terms, or clamping pivots.

This is a zero-mean additive Gaussian nuisance model independent of the state
and instrument noise. For nonlinear sensors it uses the supplied local
Jacobian, not an exact nonlinear marginalization. The caller owns units,
calibration quality, and selection of the linearization point. A registration
RMS is not a replacement for the full covariance; no calibration or scientific
validation certificate is issued here.

Run the complete synthetic example:

```sh
cargo run -p fs-assimilate --example shared_calibration
cargo test -p fs-assimilate --test factor --test robust --test factor_batch --test shared_noise
```
