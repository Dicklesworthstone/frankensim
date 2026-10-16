# Explicit noisy SensorForge acquisitions

`plan::acquisition::simulation::GaussianMeasurementSimulator` implements the
existing `MeasurementProvider`. Construct it with an explicit seed and named
`SimulationTruth` values. Only the provider owns these synthetic latent values;
`MeasurementStudy`, its alternatives and its recommendation contain no truth.
Use the ordinary `study.acquire(request, provider, cx)` and acquired-data `apply`
sequence. The legacy deterministic `Candidate.truth` example is unchanged.

Each acquired value is `truth + sqrt(requested_noise_variance) * Z` using the
existing fs-rand Philox and strict Box-Muller normal owner. It does NOT draw a
new truth from the current posterior, return a posterior prediction as data, or
ignore the positive sensor variance. These are independent additive Gaussian
instrument errors conditional on explicitly supplied fixed latent values.

The sampler hashes its version, fs-rand stream semantics, explicit seed and
complete request identity into a declared 128-bit Philox substream key. It
consumes exactly the two normal input draws at indices zero and one. The key
is a versioned hash derivation, not an injective mapping or authentication.
Re-reading one request deliberately replays one acquisition; a new request
gets a new random address. Ambient Cx randomness, provider call order and
worker scheduling do not determine the outcome. Reconstructing the same study
and simulator reproduces the closed-loop requests and acquired values.

For paired experiments the same request and seed reuse the same standardized
noise even when the supplied latent values change. Full simulator and outcome
identities still bind the latent declarations, schemas, seed and actual value.
The instrument identity is explicitly prefixed `simulated-gaussian-v1-`.

`simulator.sample(request, cx)` additionally exposes simulation-only diagnostics:
the standardized normal, exact stream key, source identity and additive noise
in the objective's difference schema. Temperature noise is a temperature delta;
the acquired value retains the absolute-temperature schema. Only the checked
outcome crosses the provider boundary. Nonfinite arithmetic, erased nonzero
perturbations and semantically inadmissible Gaussian outcomes refuse; they are
not clipped or resampled. A Gaussian model on a restricted semantic range can
therefore fail, rather than silently become a truncated distribution.

This is a reproducible synthetic instrument model, not sensor calibration,
real-device authentication, a noisy physical simulator, or statistical validation
of the decision rule. Count and compute limits retain their existing meanings.

```sh
cargo test -p fs-oed-e2e --test seeded_acquisition --test measurement
```

Regression tests exercise native-draw identity, seed and request sensitivity,
stateless replay, 2048-seed noise moments, real Kalman updates through production
acquisition, hidden-truth separation, quantity schemas and refusal/retry paths.
The moment test is a deterministic regression, not a coverage certificate.
