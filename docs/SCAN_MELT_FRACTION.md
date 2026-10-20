# Fit and constrain the material's melt state

`ScanObservable::LiquidMassFraction` observes the original constitutive chart's
liquid **mass** fraction. Values and residual scales are dimensionless fractions:
`0.8` means 80% liquid by mass. It is not a volume fraction, percentage-valued
input, or a temperature threshold. Fitting targets outside `[0,1]` refuse.

The observable works in `ScanObservation`, `ScanResponse`, and
`ScanPeakResponse`, so existing `ScanPowerStudy` L-BFGS and
`ConstrainedScanStudy` SQP can fit or constrain melt fraction while jointly
changing source powers, path coordinates, and widths. No extra forward solve
or alternate thermal model is introduced. Full and checkpointed evaluations
use the same observation chain, including intermediate-time responses and
sampled hot-spot envelopes.

For example, use the existing state-limit constructor with a row like:

```rust
ScanStateBound {
    response: ScanResponse {
        endpoint: 1,
        vertex: 0,
        observable: ScanObservable::LiquidMassFraction,
    },
    lower: Some(0.6),
    upper: Some(0.85),
    scale: 0.01,
}
```

Feasibility tolerance `tau` allows `tau * scale` physical fraction violation.
Point limits cover only their selected node and time. A fraction peak envelope
can exceed one because its positive smoothing margin is retained; it is never
clipped to manufacture feasibility. Its margin and actual sampled maximum are
reported separately, just as for temperature/enthalpy peaks.

## Executable synthetic inverse fit

```bash
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --melt-fraction
cargo run --manifest-path examples/scan_power_fit/Cargo.toml -- --melt-fraction 0.55 0.7
cargo test -p fs-conduction --test scan_fraction
cargo test --manifest-path examples/scan_power_fit/Cargo.toml
```

The unit-volume, unit-reference-density fixture starts at fraction 0.5 on a
310 K plateau, with 100 J/kg of latent enthalpy. Default targets of 0.6 after
one second and 0.8 after two seconds have an analytical schedule of 10 W then
20 W. The alternate targets imply 5 W then 15 W. The example optimizes using
the actual fraction objective, not converted enthalpy targets or hard-coded
recovered powers. It uses eight endpoints with two-step checkpoint segments,
retains the exact shutdown at two seconds, and replays the exported schedule.
A small gradient at a saturated solid/liquid state is not reported as successful
target recovery when the fraction residual remains large.

## Derivative and scope contract

Each accepted linearization retains its own node chart's `df_liquid/dh` and any
fraction-only corner/refusal metadata. Nonzero seeds at such corners reject a
classical derivative even if `T(h)` is smooth there. Unobserved zero local seeds
do not invalidate existing temperature/enthalpy derivatives. New local fraction
seeds are converted to enthalpy once; downstream history seeds are never
converted a second time. No common latent heat is assumed across materials.

This adds one derivative/refusal entry per retained mesh vertex. The existing
checkpoint/vertex allowances remain shape counts, not peak byte measurements.
Time grid, event side, geometry, constitutive charts, and reference masses stay
fixed for differentiation. Fully solid and fully liquid interiors have zero
fraction gradient, so local optimization does not guarantee escape from a
saturated plateau or a global optimum. No nucleation, hysteresis, melt flow,
unsampled extrema, or experimental calibration is claimed.

Native Rust checks could not run in the authoring session because Cargo, rustc,
DSR and RCH were absent. Independent numerical checks are equation-level
references, not evidence of native compilation or Rust optimizer execution.
