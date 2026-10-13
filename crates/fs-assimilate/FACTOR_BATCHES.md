# Factor-retaining correlated observations

`assimilate_observation_batch` now whitens a strictly positive-definite full
observation covariance and carries a single `FactorBelief` through the complete
ordered row sequence. It factors the prior once and reconstructs/revalidates the
dense posterior once at publication. The update loop is O(m n²), not O(m n³).
Observation-covariance admission/whitening and the two dense state boundaries
remain separate costs; this is not a claim of an entirely quadratic algorithm.

The exactly diagonal branch deliberately keeps `assimilate_all` unchanged,
including canonical order and bitwise replay behavior. Mismatched operator
lengths are rejected before whitening, including longer rows that the old
whitener could silently truncate.

For streaming independent or externally whitened rows, use
`factor::batch::assimilate_sequence(&factor, &rows, cx)` and retain `belief()` as
the next prior. This API uses declaration order and one cumulative work/poll
accountant for all its rows and joint residual passes. `into_parts()` moves the
factor without copying it. `FactorBelief::to_belief(cx)` provides an optional,
cancellable, checked dense publication boundary.

Correlated robust results expose `factor_diagnostics()`: scalar contraction
receipts in whitened-row order and full-batch weighted residuals evaluated at
the original and final means. The latter are NOT sums of scalar residuals
measured at intermediate means. These are numerical diagnostics of the declared
Gaussian model, not validation of that model or a general Loewner certificate.
No exact/singular observation-noise support is claimed; the existing strict
Cholesky gate still applies. The outer robust preparation/publication stages
retain their existing separate budget scopes; only the factor sequence and
checked dense conversion each have a cumulative accountant.

Tests: `cargo test -p fs-assimilate --test factor --test robust --test factor_batch`.
The new regression target compares direct block conditioning and joint
Mahalanobis residuals, the independent dense Joseph oracle, constrained prior
states, dimensional refusals, and cumulative budget exhaustion. It also
exercises a 128-state/40-observation batch beyond the dense cubic-work cap.

`groups::jacobian` constructs these full observation covariances from signed
measurement sensitivities and a multivariate nuisance model. See
`SHARED_NOISE_JACOBIANS.md` and the runnable `shared_calibration` example.

A valid correlated model can produce an exactly zero operator in a computed
whitened row, even though every physical observation operator is nonzero.
Those rows no longer invalidate the batch: their constant squared residual is
included in both joint misfits, but they do not invoke a scalar state update.
All transformed values remain available to whiten later rows. Diagnostics
expose the updated and state-independent row indices in whitened declaration
order; the physical observation audit still counts all accepted readings.
There is no small-coefficient threshold, rank certificate, or silent acceptance
of an operator erased entirely by division underflow.
