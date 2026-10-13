# SensorForge with actual measurements

`plan::measurement::MeasurementStudy` is the production ask/read/update path.
Its alternatives contain a prior, prior variance, sensor noise and acquisition
cost. There is no `truth` argument or field. The existing `Candidate` and
`run_campaign` worked-example API remains unchanged and still supplies its
caller-declared deterministic readings; it is not a noisy instrument simulator.

The study chooses which uncertain design to measure to improve the decision
about the lowest objective. It uses the existing full independent-Gaussian
expected opportunity loss (EVPI), nine-node predictive outcome quadrature, and
noise-dependent scalar Kalman model. It does not replace decision value with
variance-only ranking, a fixed fractional reduction, or a known future outcome.
Different acquired values can change the preferred design and the next request.

## Actual acquisition

Construct `MeasurementAlternative::new(name, prior, prior_variance,
sensor_noise_variance, acquisition_cost)` and then `MeasurementStudy::new(menu,
threshold, max_readings, cx)`. Priors and readings use `ObjectiveValue`; variances
are `QtyAny` in the objective's squared dimensions. The threshold uses the
objective's decision-difference schema, not an absolute-temperature schema.
Exact semantic kinds and acquisition-cost roles must agree across the menu.

`study.recommend(cx)` returns either an outcome-free `MeasurementRequest` or a
specific `MeasurementStop`. The request names the alternative, noise, cost,
estimated decision value per cost, and predicted posterior variance. A reading
limit or no-useful-action stop above threshold does not establish robustness.

Acquire the value externally, then call `request.reading(value, instrument)`
and `study.assimilate(&request, &reading, cx)`. Planning itself performs no I/O.
Requests bind the complete prior/menu/settings and current accepted history.
A request becomes stale after a successful update. A wrong-request reading,
foreign study, incompatible schema, or invalid instrument refuses. No partial
posterior, EVPI, or measurement record is published on numerical failure.
Borrowed readings remain usable for retry on the unchanged study.

`plan::acquisition::MeasurementProvider` supplies an optional real-I/O adapter:

```rust,ignore
let acquired = study.acquire(&request, &mut provider, cx)?;
// Persist the acquired data here when durable recovery is required.
let update = acquired.apply(&mut study, cx)?;
```

Acquisition does not mutate the study. `AcquisitionError::Before` means no
provider call occurred. `Provider` carries the unfulfilled request and the
provider's own error; its external side effects are not assumed reversible.
`After` carries the real acquired data as well as the refusal, including a
post-I/O cancellation/deadline failure. Reusing that data through `apply` never
calls the provider again. Providers must implement their own I/O deadlines and
cooperative cancellation, and can use the request ID as an external idempotency
key. The library does not automatically retry or refund physical acquisitions.

## Run against standard input

The example declares a labeled synthetic dimensionless model but has no hidden
truth or default outcomes. It prints the selected alternative before reading
`alternative value instrument` from standard input. A line naming a different
alternative is rejected rather than attached to the wrong request. Line storage
is bounded to 1,024 bytes. Interactive stdin is blocking; this is not a device
I/O timeout implementation.

```sh
printf 'A 10.1 synthetic-meter\nB 10.2 synthetic-meter\n' | \
  cargo run -p fs-oed-e2e --example measured_sensorforge
cargo test -p fs-oed-e2e --test measurement
```

Inspect `best_alternative`, typed `posterior_mean` / `posterior_variance`,
`updates`, `evpi`, and `history_id` between operations. Accepted records retain
actual values, instrument identities, predictions and assimilation colors.
`evpi_color()` is Estimated with infinite dispersion: no numerical remainder,
statistical model, calibration, or scientific-accuracy certificate is claimed.

## Scope

This path admits independent scalar Gaussian alternatives and independent
measurement noise. It is not a correlated-alternative model, dynamic filter,
optimal multi-step policy, noisy simulator, authenticated device protocol or
durable measurement database. The separate shared-noise stream APIs retain
calibration covariance but are not implicitly substituted into this independent
EVPI model. Menus and reading history retain the existing 256-alternative and
4,096-reading caps. Each API operation has its own bounded compute scope using
existing OED admission; this is not one lifetime compute quota or a physical
acquisition-spend cap. Provider execution is a separately controlled scope.
